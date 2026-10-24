#pragma once
#include "tnbp/ptns/function/bp.h"
#include "tnbp/ptns/function/detail/checkpoint_sha256.h"
#include <algorithm>
#include <charconv>
#include <climits>
#include "tnbp/framework/fileutility.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <type_traits>

namespace tnbp {
namespace portable_detail {
namespace fs = tnbp::file_io;
inline void require(bool b, const std::string &msg) {
  if (!b)
    throw std::runtime_error("portable checkpoint: " + msg);
}
// All ranks reach the agreement even when one rank's filesystem/codec fails.
template <class F> void phase(MPI_Comm comm, F &&f) {
  std::string error;
  try {
    f();
  } catch (const std::exception &e) {
    error = e.what();
  }
  int rank;
  MPI_Comm_rank(comm, &rank);
  int who = error.empty() ? INT_MAX : rank, first;
  MPI_Allreduce(&who, &first, 1, MPI_INT, MPI_MIN, comm);
  if (first != INT_MAX) {
    MpiBcast(error, first, comm);
    throw std::runtime_error(error);
  }
}
inline std::string file_hash(const std::string &p) {
  std::ifstream in(p, std::ios::binary);
  require(bool(in), "cannot open " + p);
  return sha256(in);
}
inline std::string text(const std::string &p) {
  require(fs::is_regular_file(p) && fs::file_size(p) <= 64 * 1024 * 1024,
          "invalid metadata file " + p);
  std::ifstream in(p, std::ios::binary);
  require(bool(in), "cannot open metadata");
  std::ostringstream out;
  out << in.rdbuf();
  require(!in.bad(), "metadata read failed");
  return out.str();
}
inline void write_text(const std::string &p, const std::string &value) {
  std::ofstream out(p, std::ios::binary | std::ios::trunc);
  out << value;
  out.close();
  require(bool(out), "write/close failed " + p);
}
inline void agree(const std::string &value, MPI_Comm comm) {
  std::string root = value;
  MpiBcast(root, 0, comm);
  phase(comm, [&] { require(root == value, "rank metadata disagreement"); });
}
inline std::uint64_t number(std::istream &in) {
  std::string s;
  require(bool(in >> s), "missing integer");
  std::uint64_t n = 0;
  auto r = std::from_chars(s.data(), s.data() + s.size(), n);
  require(r.ec == std::errc{} && r.ptr == s.data() + s.size(),
          "invalid unsigned integer");
  return n;
}
inline int site_number(std::istream &in) {
  std::string s;
  require(bool(in >> s), "missing site");
  int n = 0;
  auto r = std::from_chars(s.data(), s.data() + s.size(), n);
  require(r.ec == std::errc{} && r.ptr == s.data() + s.size(), "invalid site");
  return n;
}
inline void tag(std::istream &in, const char *want) {
  std::string s;
  require(bool(in >> s) && s == want, std::string("expected ") + want);
}
inline void end(std::istream &in) {
  in >> std::ws;
  require(in.eof(), "trailing metadata");
}
inline std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
  require(!b || a <= UINT64_MAX / b, "size overflow");
  return a * b;
}
inline std::size_t site_name_width(const std::vector<int>& sites) {
  std::size_t width=6;
  for(auto site:sites) width=std::max(width,std::to_string(site).size());
  return width;
}
inline std::size_t edge_name_width(std::size_t count) {
  return std::max<std::size_t>(6,std::to_string(count ? count-1 : 0).size());
}
inline std::string record_id(int id,std::size_t width) {
  std::ostringstream out;
  out<<std::internal<<std::setfill('0')<<std::setw(static_cast<int>(width))<<id;
  return out.str();
}
inline std::string vertex(int site,std::size_t width=6) {
  return "vertices/site-"+record_id(site,width)+".ten";
}
inline std::string message(int edge,int direction,std::size_t width=6) {
  return "messages/edge-"+record_id(edge,width)+
         (direction ? "-reverse.ten" : "-forward.ten");
}
template <class T, class = void> struct native_cpu : std::false_type {};
template <class T>
struct native_cpu<T, std::void_t<decltype(std::declval<T>().Shape())>>
    : std::true_type {};
template <class T> std::string codec() {
  using E = typename tcapi::tensor_traits<T>::elem_t;
  using R = typename tcapi::tensor_traits<T>::real_t;
  static_assert(std::is_same_v<E, R> || std::is_same_v<E, std::complex<R>>);
  static_assert(std::is_same_v<R, float> || std::is_same_v<R, double>);
  const std::uint16_t n = 1;
  bool little = *reinterpret_cast<const unsigned char *>(&n) == 1;
  std::string s =
      native_cpu<T>::value ? "gqten-native-v1" : "tcapi-cuda-record-v1";
  s += std::is_same_v<E, R> ? ":real" : ":complex";
  s += std::to_string(sizeof(R) * 8);
  if constexpr (native_cpu<T>::value)
    s += (little ? ":le:" : ":be:") + std::to_string(sizeof(std::size_t)) +
         ":" +
         std::to_string(sizeof(typename tcapi::tensor_traits<T>::bond_dim_t));
  return s;
}
struct Record {
  std::string path, hash;
  std::uint64_t length = 0;
  std::vector<std::uint64_t> shape;
  std::string line() const {
    std::ostringstream o;
    o << std::quoted(path) << ' ' << length << ' ' << hash << ' '
      << shape.size();
    for (auto d : shape)
      o << ' ' << d;
    o << '\n';
    return o.str();
  }
};
inline Record read_record(std::istream &in) {
  Record r;
  require(bool(in >> std::quoted(r.path)), "missing record path");
  r.length = number(in);
  require(bool(in >> r.hash) && r.hash.size() == 64 &&
              r.hash.find_first_not_of("0123456789abcdef") == std::string::npos,
          "invalid record hash");
  auto n = number(in);
  require(n <= 64, "tensor rank exceeds 64");
  for (std::uint64_t i = 0; i < n; ++i) {
    auto d = number(in);
    require(d > 0, "zero dimension");
    r.shape.push_back(d);
  }
  return r;
}
template <class T>
std::vector<std::uint64_t> shape(context_handle_t<T> &ctx, const T &t) {
  std::vector<std::uint64_t> s;
  for (auto d : tcapi::shape(ctx, t)) {
    require(d > 0, "invalid dimension");
    s.push_back(d);
  }
  return s;
}
// Check the native header and allocation bounds before calling the backend
// loader.
template <class T>
void check_payload(const std::string &path, const Record &r, std::uint64_t limit) {
  require(fs::is_regular_file(path) && fs::file_size(path) == r.length,
          "file length mismatch " + r.path);
  require(file_hash(path) == r.hash, "checksum mismatch " + r.path);
  using E = typename tcapi::tensor_traits<T>::elem_t;
  using D = typename tcapi::tensor_traits<T>::bond_dim_t;
  std::uint64_t count = 1;
  for (auto d : r.shape) {
    require(d <= static_cast<std::uint64_t>(std::numeric_limits<D>::max()),
            "dimension exceeds backend");
    count = mul(count, d);
  }
  auto bytes = mul(count, sizeof(E));
  require(bytes <= limit && bytes <= SIZE_MAX,
          "record allocation limit exceeded");
  std::ifstream in(path, std::ios::binary);
  auto read = [&](void *p, std::size_t n) {
    in.read(static_cast<char *>(p), n);
    require(bool(in), "truncated tensor header");
  };
  if constexpr (native_cpu<T>::value) {
    std::size_t rank;
    read(&rank, sizeof rank);
    require(rank == r.shape.size(), "native rank mismatch");
    for (auto d : r.shape) {
      D v;
      read(&v, sizeof v);
      require(v > 0 && static_cast<std::uint64_t>(v) == d,
              "native shape mismatch");
    }
    const auto header = sizeof(std::size_t) + sizeof(D) * rank + sizeof(E);
    // Scalar values live entirely in gqten's scale field.
    require(r.length == header + (rank ? bytes : 0),
            "native payload length mismatch");
  } else {
    auto uintle = [&](unsigned n) {
      std::uint64_t v = 0;
      for (unsigned i = 0; i < n; ++i) {
        unsigned char c;
        read(&c, 1);
        v |= std::uint64_t(c) << (8 * i);
      }
      return v;
    };
    char magic[8];
    read(magic, 8);
    require(std::string(magic, 8) == "TCAPIGPU" && uintle(4) == 1,
            "CUDA codec header mismatch");
    using R = typename tcapi::tensor_traits<T>::real_t;
    unsigned code = (std::is_same_v<E, R> ? 0 : 2) + (sizeof(R) == 4 ? 1 : 2);
    require(uintle(1) == code && uintle(1) == 0 && uintle(2) == 0,
            "CUDA dtype/layout mismatch");
    require(uintle(8) == r.shape.size() && uintle(8) == count,
            "CUDA rank/count mismatch");
    for (auto d : r.shape)
      require(uintle(8) == d, "CUDA shape mismatch");
    require(r.length == 32 + 8 * r.shape.size() + bytes,
            "CUDA payload length mismatch");
  }
}
} // namespace portable_detail

struct PortableCheckpoint {
  std::string directory;
  std::vector<int> sites;
  std::vector<std::pair<int, int>> edges;
  std::map<std::string, portable_detail::Record> records;
  std::string codec, state_digest;
};
namespace portable_detail {
inline std::string state_text(const PortableCheckpoint &d) {
  std::ostringstream o;
  o << "TNBP_PORTABLE_TPS 3\n"
    << std::quoted(d.codec) << '\n';
  o << "AXES incident-edge-order physical-last ket-bra\nSITES "
    << d.sites.size() << '\n';
  for (auto s : d.sites)
    o << s << '\n';
  o << "EDGES " << d.edges.size() << '\n';
  for (std::size_t e = 0; e < d.edges.size(); ++e)
    o << e << ' ' << d.edges[e].first << ' ' << d.edges[e].second << '\n';
  o << "RECORDS " << d.records.size() << '\n';
  for (auto &[key, r] : d.records)
    o << r.line();
  return o.str();
}
inline std::string manifest(const PortableCheckpoint &d) {
  std::ostringstream o;
  o << state_text(d);
  o << "END\n";
  return o.str();
}
inline void validate(const PortableCheckpoint &d) {
  const auto vw=site_name_width(d.sites),ew=edge_name_width(d.edges.size());
  require(!d.sites.empty() && std::is_sorted(d.sites.begin(), d.sites.end()) &&
              std::adjacent_find(d.sites.begin(), d.sites.end()) ==
                  d.sites.end(),
          "invalid site list");
  require(d.sites.size() <= 1000000 && d.edges.size() <= 1000000,
          "graph exceeds index range");
  std::set<std::pair<int, int>> seen;
  std::set<std::string> paths;
  for (auto s : d.sites)
    paths.insert(vertex(s,vw));
  for (std::size_t e = 0; e < d.edges.size(); ++e) {
    auto [a, b] = d.edges[e];
    require(a != b && std::binary_search(d.sites.begin(), d.sites.end(), a) &&
                std::binary_search(d.sites.begin(), d.sites.end(), b) &&
                seen.insert(std::minmax(a, b)).second,
            "invalid simple graph");
    paths.insert(message(e,0,ew));
    paths.insert(message(e,1,ew));
  }
  require(paths.size() == d.records.size(), "record count mismatch");
  for (auto &p : paths)
    require(d.records.count(p), "missing record " + p);
  for (auto s : d.sites) {
    const auto &r = d.records.at(vertex(s,vw));
    auto incident = GetSurroundingBondIndex(s, d.edges);
    require(r.shape.size() <= 64 && r.shape.size() == incident.size() + 1,
            "V axis count mismatch");
  }
  for (std::size_t e = 0; e < d.edges.size(); ++e) {
    auto [a, b] = d.edges[e];
    auto dim = [&](int s) {
      auto inc = GetSurroundingBondIndex(s, d.edges);
      auto j = std::find(inc.begin(), inc.end(), int(e)) - inc.begin();
      return d.records.at(vertex(s,vw)).shape.at(j);
    };
    auto n = dim(a);
    require(n == dim(b), "bond dimensions disagree");
    for (int direction = 0; direction < 2; ++direction)
      require(d.records.at(message(e,direction,ew)).shape ==
                  std::vector<std::uint64_t>{n, n},
              "E axis/shape mismatch");
  }
}
} // namespace portable_detail

inline PortableCheckpoint ReadTPSCheckpoint(const std::string &path) {
  using namespace portable_detail;
  PortableCheckpoint d;
  d.directory = path;
  auto m = text(fs::join(path,"manifest.txt")), c = text(fs::join(path,"COMMITTED"));
  require(c == "TNBP_COMMITTED_V1 " + sha256(m) + "\n",
          "uncommitted or corrupt manifest");
  std::istringstream in(m);
  tag(in, "TNBP_PORTABLE_TPS");
  require(number(in)==3, "unsupported schema");
  require(bool(in >> std::quoted(d.codec)), "missing codec");
  tag(in, "AXES");
  tag(in, "incident-edge-order");
  tag(in, "physical-last");
  tag(in, "ket-bra");
  tag(in, "SITES");
  auto n = number(in);
  require(n > 0 && n <= 1000000, "site count limit");
  for (std::uint64_t i = 0; i < n; ++i)
    d.sites.push_back(site_number(in));
  tag(in, "EDGES");
  n = number(in);
  require(n <= 1000000, "edge count limit");
  for (std::uint64_t i = 0; i < n; ++i) {
    require(number(in) == i, "edge ID/order mismatch");
    int a = site_number(in), b = site_number(in);
    d.edges.emplace_back(a, b);
  }
  tag(in, "RECORDS");
  n = number(in);
  require(n == d.sites.size() + 2 * d.edges.size(), "record count mismatch");
  for (std::uint64_t i = 0; i < n; ++i) {
    auto r = read_record(in);
    require(d.records.emplace(r.path, r).second, "duplicate record");
  }
  tag(in, "END");
  end(in);
  validate(d);
  d.state_digest = sha256(state_text(d));
  return d;
}
inline std::map<int, int>
PortableContiguousOwners(const std::vector<int> &sites, MPI_Comm comm) {
  int n;
  MPI_Comm_size(comm, &n);
  portable_detail::require(sites.size() >= std::size_t(n),
                           "empty ranks are unsupported");
  return ContiguousSiteOwners(sites,comm);
}
namespace portable_detail {
inline void owners_valid(const PortableCheckpoint &d,
                         const std::map<int, int> &owners, MPI_Comm comm) {
  int n;
  MPI_Comm_size(comm, &n);
  std::string signature;
  phase(comm, [&] {
    require(owners.size() == d.sites.size(), "owner coverage mismatch");
    std::vector<int> counts(n);
    for (auto s : d.sites) {
      auto it = owners.find(s);
      require(it != owners.end() && it->second >= 0 && it->second < n,
              "invalid owner");
      ++counts[it->second];
      signature += std::to_string(s) + ":" + std::to_string(it->second) + ";";
    }
    for (auto k : counts)
      require(k > 0, "empty ranks unsupported");
  });
  agree(sha256(signature), comm);
}
} // namespace portable_detail
// VerticesOnly validates the common manifest/graph and V payloads, but never
// opens, hashes or decodes an E file. The default remains a complete TPS load.
enum class PortableLoadMode { VerticesAndMessages, VerticesOnly };

template <class T> struct PortableLocalState {
  std::vector<T> V, E;
  std::vector<std::string> read_records;
  std::uint64_t payload_bytes = 0;
  std::vector<int> SiteIdx, EdgeIdx, communication_order;
  std::map<int, int> owners;
  PortableCheckpoint checkpoint;
};
template <class T>
PortableLocalState<T>
LoadTPSPortable(context_handle_t<T> &ctx, const std::string &path,
                const std::map<int, int> &owners, MPI_Comm comm,
                std::uint64_t max_record_bytes = std::uint64_t(64) * 1024 *
                                                 1024 * 1024,
                PortableLoadMode mode = PortableLoadMode::VerticesAndMessages) {
  using namespace portable_detail;
  PortableLocalState<T> state;
  int rank;
  MPI_Comm_rank(comm, &rank);
  phase(comm, [&] {
    require(mode == PortableLoadMode::VerticesAndMessages ||
                mode == PortableLoadMode::VerticesOnly,
            "invalid load mode");
    state.checkpoint = ReadTPSCheckpoint(path);
    require(state.checkpoint.codec == codec<T>(), "backend/dtype/ABI mismatch");
  });
  agree(sha256(manifest(state.checkpoint)) + ":" +
            std::to_string(static_cast<int>(mode)),
        comm);
  owners_valid(state.checkpoint, owners, comm);
  state.owners = owners;
  phase(comm, [&] {
    auto &d = state.checkpoint;
    const auto vw=site_name_width(d.sites),ew=edge_name_width(d.edges.size());
    for (auto s : d.sites)
      if (owners.at(s) == rank)
        state.SiteIdx.push_back(s);
    for (std::size_t e = 0; e < d.edges.size(); ++e)
      if (owners.at(d.edges[e].first) == rank ||
          owners.at(d.edges[e].second) == rank)
        state.EdgeIdx.push_back(int(e));
    auto load = [&](const std::string &p) {
      const auto &r = d.records.at(p);
      check_payload<T>(fs::join(path,p), r, max_record_bytes);
      require(state.payload_bytes <= UINT64_MAX - r.length,
              "read byte count overflow");
      state.payload_bytes += r.length;
      state.read_records.push_back(p);
      auto t = tcapi::load<T>(ctx, fs::join(path,p));
      require(portable_detail::shape(ctx, t) == r.shape,
              "decoded shape mismatch");
      return t;
    };
    for (auto s : state.SiteIdx)
      state.V.push_back(load(vertex(s,vw)));
    if (mode == PortableLoadMode::VerticesAndMessages) {
      for (int direction = 0; direction < 2; ++direction)
        for (auto e : state.EdgeIdx)
          state.E.push_back(load(message(e,direction,ew)));
      state.communication_order =
          MakeBPCommunicationOrder(d.edges, owners, rank);
    }
    std::clog << "PORTABLE_READ rank=" << rank << " V=" << state.V.size()
              << " E=" << state.E.size() << " bytes=" << state.payload_bytes
              << " digest=" << d.state_digest << '\n';
  });
  return state;
}

// Convenience overload without changing the existing allocation-limit argument.
template <class T>
PortableLocalState<T> LoadTPSPortable(context_handle_t<T> &ctx,
                                      const std::string &path,
                                      const std::map<int, int> &owners,
                                      MPI_Comm comm, PortableLoadMode mode) {
  return LoadTPSPortable<T>(ctx, path, owners, comm,
                            std::uint64_t(64) * 1024 * 1024 * 1024, mode);
}

// A new immutable directory is published only after every writer succeeds.
// Caller supplies a stable snapshot; this function never normalizes or repairs
// it.
template <class T>
PortableCheckpoint SaveTPSPortable(
    context_handle_t<T> &ctx, const std::string &destination,
    const std::vector<std::pair<int, int>> &edges, const std::vector<T> &V,
    const std::vector<int> &sites, const std::map<int, int> &owners,
    const std::vector<T> &E, const std::vector<int> &edge_idx,
    MPI_Comm comm) {
  using namespace portable_detail;
  PortableCheckpoint d;
  const auto target=fs::directory(destination);
  d.directory = target;
  d.edges = edges;
  d.codec = codec<T>();
  int rank, n;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &n);
  for (auto &[s, r] : owners)
    d.sites.push_back(s);
  const auto vw=site_name_width(d.sites),ew=edge_name_width(d.edges.size());
  agree(target, comm);
  agree(sha256(manifest(d)), comm);
  owners_valid(d, owners, comm);
  phase(comm, [&] {
    require(V.size() == sites.size(), "local V/site size mismatch");
    std::set<int> expected;
    for (auto &[s, r] : owners)
      if (r == rank)
        expected.insert(s);
    require(std::set<int>(sites.begin(), sites.end()) == expected &&
                sites.size() == expected.size(),
            "local site coverage mismatch");
    std::set<int> ee;
    for (std::size_t e = 0; e < edges.size(); ++e) {
      auto [a, b] = edges[e];
      require(owners.count(a) && owners.count(b), "missing edge owner");
      if (owners.at(a) == rank || owners.at(b) == rank)
        ee.insert(int(e));
    }
    require(std::set<int>(edge_idx.begin(), edge_idx.end()) == ee &&
                edge_idx.size() == ee.size() && E.size() == 2 * ee.size(),
            "local E coverage mismatch");
  });
  const std::string pending = target + ".pending";
  phase(comm, [&] {
    if (rank == 0) {
      require(!fs::exists(target), "destination already exists");
      require(fs::create_directory(pending), "pending destination exists");
      fs::create_directory(fs::join(pending,"vertices"));
      fs::create_directory(fs::join(pending,"messages"));
    }
  });
  std::string local;
  phase(comm, [&] {
    auto save = [&](const T &t, const std::string &p, bool writer) {
      Record r;
      r.path = p;
      r.shape = portable_detail::shape(ctx, t);
      // Nonwriters use a private temporary file only to check boundary
      // equality.
      auto target = writer ? fs::join(pending,p)
                           : fs::join(pending,"check-" + std::to_string(rank) + "-" + fs::filename(p));
      tcapi::save(ctx, t, target);
      r.length = fs::file_size(target);
      r.hash = file_hash(target);
      if (!writer)
        fs::remove(target);
      local += r.line();
    };
    for (std::size_t i = 0; i < V.size(); ++i)
      save(V[i], vertex(sites[i],vw), true);
    for (int direction = 0; direction < 2; ++direction)
      for (std::size_t j = 0; j < edge_idx.size(); ++j) {
        int e = edge_idx[j];
        int source = direction ? edges[e].second : edges[e].first;
        save(E[j + direction * edge_idx.size()], message(e,direction,ew),
             owners.at(source) == rank);
      }
  });
  // Only record descriptors are communicated, never Tensor payloads.
  std::map<std::string, int> copies;
  for (int r = 0; r < n; ++r) {
    std::string fragment = rank == r ? local : "";
    MpiBcast(fragment, r, comm);
    phase(comm, [&] {
      std::istringstream in(fragment);
      while (in >> std::ws && !in.eof()) {
        auto rec = read_record(in);
        auto it = d.records.find(rec.path);
        if (it != d.records.end())
          require(it->second.line() == rec.line(),
                  "boundary copies disagree " + rec.path);
        else
          d.records.emplace(rec.path, rec);
        ++copies[rec.path];
      }
    });
  }
  phase(comm, [&] {
    validate(d);
    for (auto s : d.sites)
      require(copies.at(vertex(s,vw)) == 1, "duplicate V");
    for (std::size_t e = 0; e < edges.size(); ++e)
      for (int dir = 0; dir < 2; ++dir)
        require(copies.at(message(e,dir,ew)) ==
                    (owners.at(edges[e].first) == owners.at(edges[e].second)
                         ? 1
                         : 2),
                "E copy coverage mismatch");
    if (rank == 0) {
      auto m = manifest(d);
      require(m.size() <= 64 * 1024 * 1024, "metadata size limit");
      write_text(fs::join(pending,"manifest.txt"), m);
      write_text(fs::join(pending,"COMMITTED"),
                 "TNBP_COMMITTED_V1 " + sha256(m) + "\n");
      fs::rename(pending, target);
    }
  });
  d.state_digest = sha256(state_text(d));
  return d;
}
} // namespace tnbp
