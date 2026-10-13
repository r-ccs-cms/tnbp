#include "tcapi/tcapi.h"
#include "tnbp/tnbp.h"

template <typename IntT>
void init_site_labels(const std::vector<std::pair<IntT,IntT>> & edges,
		      std::vector<IntT> & site_idx,
		      std::map<IntT,IntT> & site_to_mpi_rank,
		      MPI_Comm comm) {
  int mpi_size; MPI_Comm_size(comm,&mpi_size);
  int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
  auto sites = tnbp::GetSiteIndexFromBond(edges);
  int num_sites = static_cast<int>(sites.size());
  int site_begin = 0;
  int site_end   = num_sites;
  tnbp::get_range(mpi_size,mpi_rank,site_begin,site_end);
  site_idx.resize(site_end-site_begin);
  std::copy(sites.begin()+site_begin,
	    sites.begin()+site_end,
	    site_idx.begin());
  std::vector<int> site_to_mpi_rank_vec(num_sites,0);
  std::vector<int> site_to_mpi_rank_send(num_sites,0);
  for(const auto & i : site_idx) {
    auto it_site_adrs = std::find(sites.begin(),sites.end(),i);
    auto site_adrs = std::distance(sites.begin(),it_site_adrs);
    site_to_mpi_rank_send[site_adrs] = mpi_rank;
  }
  MPI_Allreduce(site_to_mpi_rank_send.data(),
		site_to_mpi_rank_vec.data(),
		num_sites,
		MPI_INT,
		MPI_SUM,
		comm);
  for(std::size_t i=0; i < num_sites; i++) {
    site_to_mpi_rank[sites[i]] = site_to_mpi_rank_vec[i];
  }
}

template <typename TenT, typename IntT, typename RNG>
void init_random_tensor(tcapi::context_handle_t<TenT> & ctx,
			const std::vector<std::pair<IntT,IntT>> & edges,
			std::size_t bond_dim,
			std::vector<TenT> & W,
			std::vector<IntT> & site_idx,
			std::map<IntT,IntT> & site_to_mpi_rank,
			RNG & rng,
			MPI_Comm comm) {

  using ShapeT = typename tcapi::tensor_traits<TenT>::shape_t;
  using BondDimT = typename tcapi::tensor_traits<TenT>::bond_dim_t;

  init_site_labels(edges,site_idx,site_to_mpi_rank,comm);
  W.resize(site_idx.size());
  auto itW = W.begin();
  for(auto const & i : site_idx) {
    auto bonds = tnbp::GetSurroundingBondIndex(i,edges);
    auto num_bonds = bonds.size();
    ShapeT shape_w(num_bonds,static_cast<BondDimT>(bond_dim));
    *itW = tcapi::random<TenT>(ctx,shape_w,rng);
    itW++;
  }
}

template <typename IntT, typename RNG>
void init_random_state(std::size_t p_bond_dim,
		       const std::vector<IntT> & site_idx,
		       std::vector<IntT> & site_states,
		       RNG & rng) {
  site_states.resize(site_idx.size());
  IntT dim = static_cast<IntT>(p_bond_dim);
  auto dist = std::uniform_int_distribution<IntT>(0,dim-1);
  for(auto & state : site_states) {
    state = dist(rng);
  }
}

template <typename TenT, typename IntT, typename RNG>
void init_random_tensor_state(tcapi::context_handle_t<TenT> & ctx,
	  const std::vector<std::pair<IntT,IntT>> & edges,
	  std::size_t p_bond_dim, std::size_t v_bond_dim,
	  std::vector<TenT> & W,
	  std::vector<IntT> & site_idx,
	  std::map<IntT,IntT> & site_to_mpi_rank,
	  RNG & rng, MPI_Comm comm) {
  using ShapeT = typename tcapi::tensor_traits<TenT>::shape_t;
  using BondDimT = typename tcapi::tensor_traits<TenT>::bond_dim_t;
  init_site_labels(edges,site_idx,site_to_mpi_rank,comm);
  W.resize(site_idx.size());
  auto itW = W.begin();
  for(auto const & i : site_idx) {
    auto bonds = tnbp::GetSurroundingBondIndex(i,edges);
    auto num_bonds = bonds.size();
    ShapeT shape_w(num_bonds+1,static_cast<BondDimT>(v_bond_dim));
    shape_w[num_bonds] = static_cast<BondDimT>(p_bond_dim);
    *itW++ = tcapi::random<TenT>(ctx,shape_w,rng);
  }
}

template <typename IntT>
std::vector<std::vector<IntT>> read_integer_table(const std::string& filename) {
  std::ifstream ifs(filename);
  if (!ifs) {
    throw std::runtime_error("Failed to open file: " + filename);
  }

  std::vector<std::vector<IntT>> lines;
  std::string line;

  while (std::getline(ifs, line)) {
    std::istringstream iss(line);
    std::vector<IntT> values;
    IntT x;
    while (iss >> x) {
      values.push_back(x);
    }
    lines.push_back(std::move(values));
  }
  return lines;
}

template <typename IntT>
void bcast_lines(std::vector<std::vector<IntT>> & lines,
                         int root,
                         MPI_Comm comm,
                         MPI_Datatype mpi_int_type) {

  int rank;
  MPI_Comm_rank(comm, &rank);
  std::int64_t nlines = 0;
  std::vector<std::int64_t> offsets;
  std::vector<IntT> flat;

  if (rank == root) {
    nlines = static_cast<std::int64_t>(lines.size());
    offsets.resize(nlines + 1);
    offsets[0] = 0;

    for (std::int64_t i = 0; i < nlines; ++i) {
      offsets[i + 1] = offsets[i] + static_cast<std::int64_t>(lines[i].size());
    }

    flat.reserve(offsets.back());
    for (const auto& row : lines) {
      flat.insert(flat.end(), row.begin(), row.end());
    }
  }

  MPI_Bcast(&nlines, 1, MPI_INT64_T, root, comm);

  if (rank != root) {
    offsets.resize(nlines + 1);
  }

  MPI_Bcast(offsets.data(), static_cast<int>(nlines + 1),
	    MPI_INT64_T, root, comm);

  const std::int64_t ndata = offsets.back();
  if (rank != root) {
    flat.resize(ndata);
  }

  MPI_Bcast(flat.data(), static_cast<int>(ndata),
	    mpi_int_type, root, comm);

  if (rank != root) {
    lines.resize(nlines);
    for (std::int64_t i = 0; i < nlines; ++i) {
      const auto begin = flat.begin() + offsets[i];
      const auto end   = flat.begin() + offsets[i + 1];
      lines[i].assign(begin, end);
    }
  }
}
