#include "../../apps/mcmc_sampling/typedef.h"
#include "../common/device.h"
#include "tnbp/tnbp.h"
#include <cassert>
#include <fstream>
#include <iostream>
using namespace tnbp;
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  auto comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);
  try {
    initialize_tensor_device(comm);
    if (argc < 3)
      throw std::runtime_error("mode path [output]");
    using namespace portable_detail;
    require(
        sha256("") ==
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        "SHA empty");
    require(
        sha256("abc") ==
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
        "SHA abc");
    require(
        sha256(std::string(1000000, 'a')) ==
            "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0",
        "SHA multiblock");
    ContextHandle ctx;
    tcapi::create_context(ctx);
    const std::string mode = argv[1];
    if (mode == "save" || mode == "mismatch" || mode == "wide-save") {
      const std::vector<std::pair<int, int>> edges=mode=="wide-save" ?
          std::vector<std::pair<int,int>>{{0,1},{1,9},{9,1000000}} :
          std::vector<std::pair<int,int>>{{7,20},{2,7},{-4,2}};
      const std::vector<int> all=GetSiteIndexFromBond(edges);
      auto owners = PortableContiguousOwners(all, comm);
      std::vector<Tensor> V, E;
      std::vector<int> sites, ei;
      for (auto s : all)
        if (owners.at(s) == rank) {
          sites.push_back(s);
          auto inc = GetSurroundingBondIndex(s, edges);
          tcapi::shape_t<Tensor> sh(inc.size() + 1, 2);
          V.push_back(tcapi::fill<Tensor>(ctx, sh, Elem(0.5 + 0.01 * s)));
        }
      for (std::size_t i = 0; i < edges.size(); ++i)
        if (owners.at(edges[i].first) == rank ||
            owners.at(edges[i].second) == rank)
          ei.push_back(i);
      for (int dir = 0; dir < 2; ++dir)
        for (int e : ei) {
          auto x = 0.25 + e + dir * 0.75;
          if (mode == "mismatch")
            x += rank;
          E.push_back(
              tcapi::fill<Tensor>(ctx, tcapi::shape_t<Tensor>{2, 2}, Elem(x)));
        }
      auto d = SaveTPSPortable(ctx, argv[2], edges, V, sites, owners, E, ei,
                               comm);
      if (rank == 0)
        std::cout << "DIGEST " << d.state_digest << '\n';
    } else {
      PortableCheckpoint d;
      phase(comm, [&] { d = ReadTPSCheckpoint(argv[2]); });
      const auto vw=site_name_width(d.sites),ew=edge_name_width(d.edges.size());
      auto owners = PortableContiguousOwners(d.sites, comm);
      if (mode == "reverse") {
        int size;
        MPI_Comm_size(comm, &size);
        for (auto &[s, r] : owners)
          r = size - 1 - r;
      }
      if (mode == "bad-owner")
        owners.erase(owners.begin());
      auto state = LoadTPSPortable<Tensor>(
          ctx, argv[2], owners, comm, mode == "limit" ? 1 : UINT64_MAX,
          mode == "vertices-only" ? PortableLoadMode::VerticesOnly
                                  : PortableLoadMode::VerticesAndMessages);
      if (mode == "vertices-only") {
        require(state.E.empty() && state.communication_order.empty(),
                "E constructed in vertices-only mode");
        require(state.read_records.size() == state.V.size(),
                "unexpected payload reads");
        std::uint64_t bytes = 0;
        for (auto site : state.SiteIdx)
          bytes += d.records.at(vertex(site,vw)).length;
        require(state.payload_bytes == bytes, "unexpected payload bytes");
        for (const auto &path : state.read_records)
          require(path.rfind("vertices/", 0) == 0,
                  "E file read in vertices-only mode");
      }
      // Compare bytes from actual decoded tensors, not only manifest hashes.
      for (std::size_t i = 0; i < state.V.size(); ++i) {
        std::ostringstream o(std::ios::binary);
        tcapi::save(ctx, state.V[i], o);
        require(sha256(o.str()) == d.records.at(vertex(state.SiteIdx[i],vw)).hash,
                "V round trip");
      }
      for (int dir = 0; mode != "vertices-only" && dir < 2; ++dir)
        for (std::size_t j = 0; j < state.EdgeIdx.size(); ++j) {
          std::ostringstream o(std::ios::binary);
          tcapi::save(ctx, state.E[j + dir * state.EdgeIdx.size()], o);
          require(sha256(o.str()) ==
                      d.records.at(message(state.EdgeIdx[j],dir,ew)).hash,
                  "directed E round trip");
        }
      if (argc > 3) {
        auto copy = SaveTPSPortable(ctx, argv[3], d.edges, state.V,
                                    state.SiteIdx, owners, state.E,
                                    state.EdgeIdx, comm);
        require(copy.state_digest==d.state_digest,"resave digest");
      }
      if (rank == 0)
        std::cout << "DIGEST " << d.state_digest << " VERIFIED_BYTES\n";
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    MPI_Abort(comm, 1);
  }
  MPI_Finalize();
  return 0;
}
