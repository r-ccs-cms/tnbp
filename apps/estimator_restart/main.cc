// main.cc
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "qasm/any.h"
#include "qasm/utility.h"
#include "pauli/sparse_pauli.h"
#include "tnbp/tnbp.h"

#include "typedef.h"
#include "../cuda_device.h"
#include "tnbp/framework/timestamp.h"
#include "option.h"
#include "checkpoint.h"
#include "../site_ownership.h"

int main(int argc, char * argv[]) {

  int mpi_err = MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int mpi_master = 0;
  int mpi_size; MPI_Comm_size(comm,&mpi_size);
  int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
  
  try {
  initialize_tensor_device(comm);

  Option options =generate_options(argc,argv);
  if(mpi_rank==0) cout_options(options);

  if(!options.portable_save.empty() && !options.loadname.empty())
    throw std::invalid_argument("legacy import requires a separate converter");
  if(!options.portable_save.empty() && options.portable_load.empty() && options.step_start!=0)
    throw std::invalid_argument("fresh portable evolution must start at step zero");
  std::string qasm_source;
  if( mpi_rank == mpi_master ) {
    qasm_source = qasm::read_all_file(options.circuit);
  }
  tnbp::MpiBcast(qasm_source,mpi_master,comm);
  qasm::Program qasm_program = qasm::parse_any(qasm_source);

  ContextHandle ctx;
  tcapi::create_context(ctx);
  std::vector<std::pair<int,int>> edges;
  
  if( options.backend == std::string("ibm_kobe") ) {
    edges = tnbp::bond_ibm_kobe();
  } else if ( options.backend == std::string("default") ) {
    edges = tnbp::EdgesFromQasm(qasm_program);
  }

  // The saved tensor axes define the graph, independently of the new circuit.
  tnbp::PortableCheckpoint descriptor;
  if(!options.portable_load.empty()) {
    tnbp::portable_detail::phase(comm,[&] {
      descriptor=tnbp::ReadTPSCheckpoint(options.portable_load);
      if(options.backend!="default" && descriptor.edges!=edges)
        throw std::runtime_error("portable graph differs from requested backend graph/order");
      edges=descriptor.edges;
      for(const auto& instruction:qasm_program.instructions)
        for(const auto& q:instruction.qubits)
          if(std::find(descriptor.sites.begin(),descriptor.sites.end(),q.index)==descriptor.sites.end())
            throw std::runtime_error("circuit references site outside checkpoint graph");
    });
  }

  if( mpi_rank == mpi_master ) {
    std::cout << " " << tnbp::make_timestamp()
	      << " edges:";
    for(const auto & [u,v] : edges) {
      std::cout << " (" << u << "," << v << ")";
    }
    std::cout << std::endl;
    
    std::cout << " " << tnbp::make_timestamp()
	      << " Start TPO construction for quantum circuit " << std::endl;
  }

  std::vector<std::vector<Tensor>> TPO;
  TPO = tnbp::QasmToTPO<Tensor>(ctx,qasm_program,edges);
  if( mpi_rank == mpi_master ) {
    std::cout << " " << tnbp::make_timestamp()
	      << " barrier-based layering of QasmToTPO was performed " << std::endl;
  }

  std::vector<int> qubit = tnbp::GetSiteIndexFromBond(edges);
  std::vector<Tensor> V;
  std::vector<int> SiteIdx;
  std::map<int,int> Site_To_MpiRank;
  std::vector<Tensor> E;
  std::vector<int> EdgeIdx;
  std::vector<int> pdim(qubit.size(),2);
  std::vector<Tensor> F;

  tnbp::BPCheckpointMetadata bp_metadata;
  bp_metadata.source_algorithm="fresh";
  if(!options.portable_load.empty() && !options.loadname.empty())
    throw std::invalid_argument("choose portable_load or loadname");
  auto estimator_settings = [&] {
    std::ostringstream out;out<<std::setprecision(17)<<options.max_bond_dim<<' '
      <<options.max_bp_iterations<<' '<<options.bp_tolerance<<' '<<options.sv_min<<' '<<options.truncation_error;
    return out.str();
  };
  std::size_t portable_next_step=0;
  if(!options.portable_load.empty()) {
    tnbp::portable_detail::phase(comm,[&] {
      auto resume=estimator_checkpoint::read(options.portable_load,descriptor.state_digest);
      bool same_circuit=false;
      if(resume) {
        same_circuit=resume->circuit_sha256==tnbp::portable_detail::sha256(qasm_source);
        if(same_circuit && resume->settings!=estimator_settings())
          throw std::runtime_error("estimator resume conditions differ");
        portable_next_step=same_circuit ? resume->next_step : 0;
      }
      if(options.step_start_set && options.step_start!=portable_next_step)
        throw std::runtime_error("step_start conflicts with inferred start step");
      if(mpi_rank==0)
        std::cout<<"PORTABLE_RESTART circuit="<<(resume ? (same_circuit ? "same" : "changed") : "unspecified")
                 <<" step_start="<<portable_next_step<<std::endl;
    });
    auto owners=app_ownership::select(descriptor.sites,options.ownerfile,comm,true);
    auto loaded=tnbp::LoadTPSPortable<Tensor>(ctx,options.portable_load,owners,comm);
    V=std::move(loaded.V); E=std::move(loaded.E); SiteIdx=std::move(loaded.SiteIdx);
    EdgeIdx=std::move(loaded.EdgeIdx); Site_To_MpiRank=std::move(loaded.owners);
  } else if( !options.loadname.empty() ) {
    std::string filename = tnbp::make_rank_filename(options.loadname,mpi_rank);
    std::ifstream ifs(filename,std::ios::binary);
    std::cout << " " << tnbp::make_timestamp()
	      << " Start loading of TPS data at rank " << mpi_rank << std::endl;
    tnbp::LoadTPS(ctx,ifs,V,SiteIdx,Site_To_MpiRank,E,EdgeIdx,&bp_metadata);
    bp_metadata.source_checkpoint=filename;
    std::cout << " " << tnbp::make_timestamp()
	      << " Load TPS data at rank " << mpi_rank << std::endl;
  } else {
    Site_To_MpiRank=app_ownership::select(qubit,options.ownerfile,comm);
    tnbp::InitTensorProductState(ctx,edges,pdim,
				 V,SiteIdx,Site_To_MpiRank,
				 E,EdgeIdx,comm);
  }
  

  if(!options.loadname.empty() && !options.ownerfile.empty()) {
    auto requested=app_ownership::select(qubit,options.ownerfile,comm);
    app_ownership::check_legacy(Site_To_MpiRank,requested,comm);
  }
  app_ownership::print(Site_To_MpiRank,options.ownerfile.empty() ?
    (options.loadname.empty() ? "automatic" : "checkpoint") : options.ownerfile,comm);

  Real tolerance;
  std::vector<BondDim> res_bond_dim;
  std::vector<Real> res_truncation_error;

  if( mpi_rank == mpi_master ) {
    std::cout << " " << tnbp::make_timestamp()
	      << " Start circuit contractions " << std::endl;
  }

  /**
     Setup for measurements
  */
  auto spo = pauli::load_sparse_pauli_op<Elem>(options.sparsepauli);
  std::vector<Tensor> exOp;
  std::vector<std::vector<int>> exSite;
  tnbp::SparsePauliToTensorOp(ctx,spo,exOp,exSite);
  if( mpi_rank == mpi_master ) {
    std::cout << " " << tnbp::make_timestamp()
	      << " number of measrement operators = " << exOp.size() << std::endl;
  }
  std::vector<Tensor> exOp_one;
  std::vector<int> exSite_one;
  std::vector<Tensor> exOp_two;
  std::vector<std::pair<int,int>> exPair_two;
  tnbp::ClassifyOps(ctx,edges,exOp,exSite,
		    exOp_one,exSite_one,
		    exOp_two,exPair_two);
  if( mpi_rank == mpi_master ) {
    std::cout << " " << tnbp::make_timestamp()
	      << " generated single site operator " << std::endl;
  }

  size_t step_start = 0;
  size_t step_end = 1;
  size_t barrier_start = 0;
  size_t barrier_end = TPO.size();
  if( options.step_start != 0 ) {
    step_start = options.step_start;
  }
  if( options.step_end != 0 ) {
    step_end = options.step_end;
  }

  if(!options.portable_load.empty()) step_start=portable_next_step;
  if((!options.portable_save.empty()||!options.portable_load.empty()) && step_end<step_start)
    throw std::invalid_argument("step_end precedes checkpoint next step");
  const auto bp_order=tnbp::MakeBPCommunicationOrder(edges,Site_To_MpiRank,mpi_rank);
  for(size_t t=step_start; t < step_end; t++) {
    for(size_t m=barrier_start; m < barrier_end; m++) {
      /**
	 Attach operator tensor to site tensor
      */
      tnbp::AbsorbTPO(
		      ctx,TPO[m],edges,
		      V,SiteIdx,Site_To_MpiRank,E,EdgeIdx,
		      comm);
      
      /**
	 Perform belief propagation loop
      */
      for(size_t k=0; k < options.max_bp_iterations; k++) {
      
	tnbp::UpdateBPMessages(ctx,edges,V,SiteIdx,E,EdgeIdx,comm,F);
	tnbp::ExchangeBPMessages(ctx,edges,Site_To_MpiRank,EdgeIdx,bp_order,comm,F);
	// Compare before committing the complete synchronous update.
	tnbp::BeliefPropagationDifference(ctx,edges,Site_To_MpiRank,
                                            E,F,EdgeIdx,comm,tolerance);
      E.swap(F);
	if( mpi_rank == mpi_master ) {
	  std::cout << " " << tnbp::make_timestamp()
		    << " error from belief propagation condition at step "
		    << t << " and barrier " << m << " = "
		    << tolerance << std::endl;
	}
	if( tolerance < options.bp_tolerance ) {
	  break;
	}
      }
    
      /**
	 Perform truncation
      */
      tnbp::Truncation(ctx,edges,
		       V,SiteIdx,Site_To_MpiRank,
		       E,EdgeIdx,comm,
		       options.max_bond_dim,
		       options.sv_min,
		       options.truncation_error,
		       res_bond_dim,
		       res_truncation_error);
      for(size_t k=0; k < EdgeIdx.size(); k++) {
	std::cout << " " << tnbp::make_timestamp()
		  << " truncation error for edge ("
		  << edges[EdgeIdx[k]].first << ","
		  << edges[EdgeIdx[k]].second << ") = "
		  << res_truncation_error[k] << " bond dim = "
		  << res_bond_dim[k] << std::endl;
      }
      
      /**
	 Measurement step
      */
      auto it_meas = std::find(options.measurement_barrier.begin(),
			       options.measurement_barrier.end(),
			       m);
      if( it_meas != options.measurement_barrier.end() ) {
	
	if( mpi_rank == mpi_master ) {
	  std::cout << " " << tnbp::make_timestamp()
		    << " start measurements " << std::endl;
	}
	
	std::vector<Elem> exVal_one =
	  tnbp::Measure(ctx,edges,V,SiteIdx,Site_To_MpiRank,
			E,EdgeIdx,comm,exSite_one,exOp_one);
	
	std::cout.precision(16);
	for(int site_address=0; site_address < exSite_one.size(); site_address++) {
	  if( mpi_rank == Site_To_MpiRank[exSite_one[site_address]] ) {
	    std::cout << " Expectation value of single site operator for step "
		      << t << " barrier " << m << " at site "
		      << exSite_one[site_address]
		      << " " << GetReal(exVal_one[site_address]) << std::endl;
	  }
	}
	
	std::vector<Elem> exVal_two_rank =
	  tnbp::Measure(ctx,edges,V,SiteIdx,Site_To_MpiRank,
			E,EdgeIdx,comm,exPair_two,exOp_two);
	std::vector<Elem> exVal_two(exPair_two.size());
	
	MPI_Datatype DataT = tnbp::GetMpiType<Elem>();
	MPI_Allreduce(exVal_two_rank.data(),exVal_two.data(),
		      exVal_two.size(),DataT,MPI_SUM,comm);
	
	if( mpi_rank == mpi_master ) {
	  for(int edge_address=0; edge_address < exPair_two.size(); edge_address++) {
	    std::cout << " Expectation value of two site operator for step "
		      << t << " barrier " << m << " at sites ("
		      << exPair_two[edge_address].first
		      << "," << exPair_two[edge_address].second
		      << ") "
		      << GetReal(exVal_two[edge_address]) << std::endl;
	  }
	}
      } // end measurements
    }
  }

  if( !options.savename.empty() ) {
    std::string filename = tnbp::make_rank_filename(options.savename,mpi_rank);
    std::ofstream ofs(filename,std::ios::binary);
    tnbp::SaveTPS(ctx,ofs,V,SiteIdx,Site_To_MpiRank,E,EdgeIdx,bp_metadata);
  }

  if(!options.portable_save.empty()) {
    auto checkpoint=estimator_checkpoint::save(options.portable_save,
      {"",tnbp::portable_detail::sha256(qasm_source),estimator_settings(),step_end},comm,
      [&](const std::string& state) {
        return tnbp::SaveTPSPortable(ctx,state,edges,V,SiteIdx,Site_To_MpiRank,E,EdgeIdx,comm);
      });
    if(mpi_rank==0) std::cout<<"PORTABLE_SAVED digest="<<checkpoint.state_digest<<std::endl;
  }
  } catch(const std::exception& error) {
    std::cerr<<"estimator: "<<error.what()<<std::endl;MPI_Abort(comm,1);
  }
  MPI_Finalize();
  
  return 0;
}
