#include "tci/tci.h"
#include "tnbp/tnbp.h"

template <typename TenT, typename IntT>
void init_random_tensor(tci::context_handle_t<TenT> & ctx,
			const std::vector<std::pair<IntT,IntT>> & edges,
			size_t bond_dim,
			std::vector<TenT> & W,
			std::vector<IntT> & site_idx,
			std::map<IntT,IntT> & site_to_mpi_rank,
			uint32_t seed,
			MPI_Comm comm) {

  using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
  using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;

  std::mt19937 rng(seed);
  
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
  W.resize(site_end-site_begin);
  auto itW = W.begin();
  for(auto const & i : site_idx) {
    auto bonds = tnbp::GetSurroundingBondIndex(i,edges);
    auto num_bonds = bonds.size();
    ShapeT shape_w(num_bonds,static_cast<BondDimT>(bond_dim));
    *itW = tci::random<TenT>(ctx,shape_w,rng);
    itW++;
  }
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
  for(size_t i=0; i < num_sites; i++) {
    site_to_mpi_rank[sites[i]] = site_to_mpi_rank_vec[i];
  }
}
			   
