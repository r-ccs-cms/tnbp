#include <vector>
#include <utility>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <mpi.h>

template <class IntT>
std::vector<std::pair<IntT, IntT>>
read_edges_from_file(const std::string& filename, bool deduplicate) {
  std::ifstream ifs(filename.c_str());
  if (!ifs) {
    throw std::runtime_error("failed to open file: " + filename);
  }
  std::vector<std::pair<IntT, IntT>> edges;
  std::string line;
  int line_no = 0;
  while (std::getline(ifs, line)) {
    ++line_no;
    // delete line with '#'
    std::string::size_type pos = line.find('#');
    if (pos != std::string::npos) {
      line.erase(pos);
    }
    std::istringstream iss(line);
    IntT u, v;
    if (!(iss >> u)) {
      continue;
    }
    if (!(iss >> v)) {
      throw std::runtime_error(
	    "parse error at line " + std::to_string(line_no) +
	    ": expected two values");
    }
    // check whether non-meaningful token
    std::string extra;
    if (iss >> extra) {
      throw std::runtime_error(
	    "parse error at line " + std::to_string(line_no) +
	    ": too many tokens");
    }
    if (u == v) {
      throw std::runtime_error(
	    "invalid self-loop at line " + std::to_string(line_no));
    }
    edges.push_back(tnbp::make_edge(u, v));
  }

  if( deduplicate ) {
    std::sort(edges.begin(),edges.end());
    edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
  }
  
  return edges;
}

template <class IntT>
void bcast_edges(std::vector<std::pair<IntT, IntT>>& edges,
                 int root,
                 MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  
  // distribute size
  unsigned long long num_edges_ull = 0;
  if (rank == root) {
    num_edges_ull = static_cast<unsigned long long>(edges.size());
  }
  
  MPI_Bcast(&num_edges_ull, 1, MPI_UNSIGNED_LONG_LONG, root, comm);

  // count in MPI_Bcast is int: check size
  if (num_edges_ull >
      static_cast<unsigned long long>(std::numeric_limits<int>::max() / 2)) {
    throw std::runtime_error("too many edges for single MPI_Bcast");
  }

  const std::size_t num_edges = static_cast<std::size_t>(num_edges_ull);
  const int buf_count = static_cast<int>(2 * num_edges);
  
  std::vector<long long> buf;
  buf.resize(static_cast<std::size_t>(buf_count));
  
  if (rank == root) {
    for (std::size_t k = 0; k < num_edges; ++k) {
      buf[2 * k    ] = static_cast<long long>(edges[k].first);
      buf[2 * k + 1] = static_cast<long long>(edges[k].second);
    }
  }

  MPI_Bcast(buf.data(), buf_count, MPI_LONG_LONG, root, comm);
  
  if (rank != root) {
    edges.resize(num_edges);
    for (std::size_t k = 0; k < num_edges; ++k) {
      edges[k].first  = static_cast<IntT>(buf[2 * k    ]);
      edges[k].second = static_cast<IntT>(buf[2 * k + 1]);
    }
  }
}


template <typename IntT>
std::string step_to_string(IntT step, int width = 6) {
  static_assert(std::is_integral<IntT>::value, "IntT must be an integral type.");
  std::ostringstream oss;
  oss << std::setw(width) << std::setfill('0') << step;
  return oss.str();
}
  
template <typename StepT>
std::string make_step_filename(const std::string & prefix, StepT step,
			       const std::string & ext = ".txt",
			       int width = 6) {
  return prefix + step_to_string(step,width) + ext;
}

template <typename IntT, typename StepT>
void write_edges_to_file(const std::string & filename,
			 const std::vector<std::pair<IntT,IntT>> & edges,
			 StepT step) {

  std::ofstream ofs(filename.c_str());
  ofs << "# graph at step " << step << "\n";
  for(auto const & [x,y] : edges) {
    ofs << x << " " << y << "\n";
  }
  ofs.close();
}
