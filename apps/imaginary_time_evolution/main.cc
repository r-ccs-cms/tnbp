#include "tnbp/tnbp.h"
#include "tnbp/framework/graphio.h"
#include "typedef.h"
#include "checkpoint.h"
#include "../cuda_device.h"
#include "../site_ownership.h"
#include "../mcmc_sampling/initial_states.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>

namespace ite {
template<class... Args> void log(bool input,const Args&... args) {
  std::ostringstream out;out<<std::setprecision(17);
  if(input)out<<"# ";out<<tnbp::make_timestamp()<<' ';(out<<...<<args);std::cout<<out.str()<<std::endl;
}
template<class F> std::string root_text(MPI_Comm comm,F f) {
  int rank;MPI_Comm_rank(comm,&rank);std::string text,error;
  if(rank==0)try{text=f();}catch(const std::exception& e){error=e.what();}
  tnbp::MpiBcast(error,0,comm);if(!error.empty())throw std::runtime_error(error);
  tnbp::MpiBcast(text,0,comm);return text;
}
struct Measurements {
  std::vector<Tensor> one,two;
  std::vector<int> sites;
  std::vector<std::pair<int,int>> pairs;
  std::vector<double> c_one,c_two;
  double identity=0;
  std::size_t unsupported=0;
};
inline Measurements measurement_ops(ContextHandle& ctx,const tnbp::PauliEvolutionPlan& p) {
  Measurements m;m.identity=p.hamiltonian.identity_coefficient;
  for(const auto& term:p.hamiltonian.terms) {
    const auto& ops=term.local_ops;
    bool supported=ops.size()==1;
    if(ops.size()==2) supported=tnbp::contains_edge(p.graph.edges,static_cast<int>(ops[0].qubit),static_cast<int>(ops[1].qubit));
    if(!supported){++m.unsupported;continue;}
    std::vector<Tensor> tensors;std::vector<std::vector<int>> sites;
    tnbp::SparsePauliToTensorOp<Tensor>(ctx,{{Elem(1),term.label}},tensors,sites);
    if(ops.size()==1) {m.one.push_back(std::move(tensors[0]));m.sites.push_back(sites[0][0]);m.c_one.push_back(term.coefficient);}
    else {m.two.push_back(std::move(tensors[0]));m.pairs.push_back({sites[0][0],sites[0][1]});m.c_two.push_back(term.coefficient);}
  }
  return m;
}
struct Window {
  std::size_t start;
  double error_min=std::numeric_limits<double>::infinity(),error_max=-std::numeric_limits<double>::infinity();
  void add(double e) {require(std::isfinite(e),"nonfinite truncation error");error_min=std::min(error_min,e);error_max=std::max(error_max,e);}
};
inline std::string number(double x) {
  if(!std::isfinite(x))return "N/A";std::ostringstream s;s<<std::setprecision(17)<<x;return s.str();
}
inline void summary(ContextHandle& ctx,const Config& cfg,const tnbp::PauliEvolutionPlan& p,
    const Measurements& m,const std::vector<Tensor>& V,const std::vector<int>& sites,
    const std::map<int,int>& owners,const std::vector<Tensor>& E,const std::vector<int>& edge_idx,
    Window& window,std::size_t next,MPI_Comm comm) {
  int rank;MPI_Comm_rank(comm,&rank);
  auto one=tnbp::Measure(ctx,p.graph.edges,V,sites,owners,E,edge_idx,comm,m.sites,m.one);
  auto two=tnbp::Measure(ctx,p.graph.edges,V,sites,owners,E,edge_idx,comm,m.pairs,m.two);
  Elem local{};
  for(std::size_t i=0;i<one.size();++i)local+=m.c_one[i]*one[i];
  for(std::size_t i=0;i<two.size();++i)local+=m.c_two[i]*two[i];
  Elem total{};MPI_Reduce(&local,&total,1,tnbp::GetMpiType<Elem>(),MPI_SUM,0,comm);
  double dmin=std::numeric_limits<double>::infinity(),dmax=-std::numeric_limits<double>::infinity();
  for(std::size_t j=0;j<sites.size();++j) {
    const auto shape=tcapi::shape(ctx,V[j]);const auto& incident=p.graph.incident_edge_ids[sites[j]];
    for(std::size_t k=0;k<incident.size();++k)if(p.graph.edges[incident[k]].first==sites[j]) {
      dmin=std::min(dmin,double(shape[k]));dmax=std::max(dmax,double(shape[k]));
    }
  }
  double mins[2]={dmin,window.error_min},maxs[2]={dmax,window.error_max},global_min[2],global_max[2];
  MPI_Reduce(mins,global_min,2,MPI_DOUBLE,MPI_MIN,0,comm);MPI_Reduce(maxs,global_max,2,MPI_DOUBLE,MPI_MAX,0,comm);
  if(rank==0) {
    total+=m.identity;require(std::isfinite(total.real()) && std::isfinite(total.imag()),"nonfinite measured energy");
    log(false,"SUMMARY step=",next," tau=",elapsed_time(cfg.schedule,next)," time_step=",next ? number(step_size(cfg.schedule,next-1)):"N/A",
      m.unsupported ? " energy=N/A partial_energy=" : " energy=",total.real()," energy_imag=",total.imag()," unsupported_terms=",m.unsupported,
      " bond_dim_min=",number(global_min[0])," bond_dim_max=",number(global_max[0]),
      " trunc_error_min=",number(global_min[1])," trunc_error_max=",number(global_max[1]),
      " window_start=",window.start," window_end=",next);
  }
  window=Window{next};
}
inline void help() {
  const std::vector<std::string> lines={
    "imaginary_time_evolution: BP imaginary-time evolution of a Pauli Hamiltonian",
    "--inputedge FILE --sparse_pauli FILE: ordered graph and real [imag] PauliString rows",
    "--initial_state zero|plus (default zero), or --initial_statefile FILE (site bit rows)",
    "--time_step DT --num_steps N, or --schedule FILE (DT N rows)",
    "--trotter_order 1|2 (default 2); one step is the complete Trotter sequence",
    "--max_bp_iterations 50 --bp_tolerance 1e-8 --max_bond_dim 100 --sv_min 1e-8 --truncation_error 1e-8",
    "--measure_every K (default 1), --step_end K (exclusive absolute endpoint)",
    "--portable_load DIR --portable_save NEW_DIR --restart_mode auto|initial (default auto)",
    "--ownerfile FILE (site rank rows); omitted uses balanced contiguous ownership",
    "Resume restores problem/schedule/settings; explicit changes require --restart_mode initial",
    "Measurement supports one-site and graph-edge two-site terms; no extra measurement BP",
    "Normal termination always measures; periodic checkpoints and mid-step resume are unsupported"};
  for(const auto& line:lines)log(true,line);
}
} // namespace ite

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);MPI_Comm comm=MPI_COMM_WORLD;int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
  try {
    const auto o=ite::options(argc,argv);
    if(o.help){if(rank==0)ite::help();MPI_Finalize();return 0;}
    initialize_tensor_device(comm);
    std::vector<std::pair<int,int>> edges;
    tnbp::PortableCheckpoint descriptor;
    if(o.has("portable_load")) {
      tnbp::portable_detail::phase(comm,[&]{descriptor=tnbp::ReadTPSCheckpoint(o.get("portable_load"));});edges=descriptor.edges;
    }
    std::size_t next=0;double operator_log_scale=0;bool continuing=false;
    auto encoded=ite::root_text(comm,[&] {
      std::optional<ite::Resume> saved;
      if(o.has("portable_load")) {
        saved=ite::read_resume(o.get("portable_load"),descriptor.state_digest);
        if(saved) {
          auto old_plan=tnbp::MakePauliEvolutionPlan(saved->config.h,edges);
          ite::require(old_plan.plan_digest==saved->plan_digest,"saved plan digest/version mismatch");
        }
      }
      continuing=bool(saved) && o.get("restart_mode","auto")!="initial";
      ite::Config cfg;
      if(continuing) {
        ite::require(saved->operator_weights==tnbp::PauliExponentialWeights,
          "saved operator weights differ; use --restart_mode initial for a new calculation");
        cfg=saved->config;auto specified=cfg;
        if(o.has("sparse_pauli"))specified.h=ite::read_hamiltonian(o.get("sparse_pauli"));
        if(o.has("schedule") || o.has("time_step"))specified.schedule=ite::read_schedule(o);
        specified.settings=ite::settings(o,cfg.settings);
        ite::require(ite::encode_config(cfg)==ite::encode_config(specified),"resume conditions differ; use --restart_mode initial for a new calculation");
        next=saved->next_step;operator_log_scale=saved->operator_log_scale;
      } else {
        ite::require(o.has("sparse_pauli"),"new calculation requires sparse_pauli");
        cfg.h=ite::read_hamiltonian(o.get("sparse_pauli"));cfg.schedule=ite::read_schedule(o);cfg.settings=ite::settings(o);
      }
      std::ostringstream out;out<<std::setprecision(17)<<next<<' '<<operator_log_scale<<' '<<continuing<<' '<<std::quoted(ite::encode_config(cfg));return out.str();
    });
    std::istringstream envelope(encoded);std::string config_text;envelope>>next>>operator_log_scale>>continuing>>std::quoted(config_text);
    auto cfg=ite::decode_config(config_text);
    const auto measure_every=ite::count(o.get("measure_every","1"));ite::require(measure_every>0,"measure_every must be positive");
    const auto end=o.has("step_end") ? ite::count(o.get("step_end")):ite::total_steps(cfg.schedule);
    ite::require(next<=end && end<=ite::total_steps(cfg.schedule),"step_end outside saved step/schedule range");
    // Only metadata is shared here; no synchronization is introduced per layer.
    if(!o.has("portable_load")) {
      ite::root_text(comm,[&]{ite::require(o.has("inputedge"),"fresh calculation requires inputedge");edges=tnbp::read_edges_from_file<int>(o.get("inputedge"),false);return std::string();});
      tnbp::bcast_edges(edges,0,comm);
    } else if(o.has("inputedge")) {
      ite::root_text(comm,[&]{ite::require(tnbp::read_edges_from_file<int>(o.get("inputedge"),false)==edges,"graph/order differs from portable state");return std::string();});
    }
    auto plan=tnbp::MakePauliEvolutionPlan(cfg.h,edges);const auto& all_sites=plan.graph.sites;
    if(o.has("portable_load"))ite::require(all_sites==descriptor.sites,"checkpoint site IDs must be 0..N-1");
    auto owners=app_ownership::select(all_sites,o.get("ownerfile"),comm,o.has("portable_load"));
    if(rank==0) {
      ite::log(true,"app: imaginary_time_evolution");for(const auto& [key,value]:o.values)ite::log(true,key,": ",value);
      ite::log(true,"initialization: ",o.has("portable_load") ? "checkpoint" : o.has("initial_statefile") ? "basis-file" : o.get("initial_state","zero"));
      ite::log(true,"site_ownership: ",o.get("ownerfile","automatic"));
      ite::log(true,"operator_weights: ",tnbp::PauliExponentialWeights);
      ite::log(true,"config: ",ite::encode_config(cfg));ite::log(true,"plan_digest: ",plan.plan_digest," layers: ",plan.layers.size());
      ite::log(true,"step_start: ",next," step_end: ",end," continuing: ",continuing," measure_every: ",measure_every);
      for(std::size_t e=0;e<edges.size();++e)ite::log(true,"edge ",e,": ",edges[e].first," ",edges[e].second);
      for(int r=0;r<size;++r) {std::ostringstream row;for(auto [site,owner]:owners)if(owner==r)row<<' '<<site;ite::log(true,"MPI rank ",r,":",row.str());}
    }
    ContextHandle ctx;tcapi::create_context(ctx);
    std::vector<Tensor> V,E,F;std::vector<int> sites,edge_idx;
    if(o.has("portable_load")) {
      auto loaded=tnbp::LoadTPSPortable<Tensor>(ctx,o.get("portable_load"),owners,comm);
      V=std::move(loaded.V);E=std::move(loaded.E);sites=std::move(loaded.SiteIdx);edge_idx=std::move(loaded.EdgeIdx);
      for(std::size_t i=0;i<V.size();++i)ite::require(tcapi::shape(ctx,V[i]).back()==2,"checkpoint physical dimension must be two");
    } else {
      tnbp::InitTensorProductState(ctx,all_sites,edges,std::vector<int>(all_sites.size(),2),V,sites,owners,E,edge_idx,comm);
      std::vector<std::vector<int>> bits;
      if(o.has("initial_statefile")) {
        auto text=ite::root_text(comm,[&]{bits=read_initial_states(o.get("initial_statefile"),all_sites,1);std::ostringstream out;for(const auto& row:bits){ite::require(row[0]<=1,"initial_statefile bit must be 0 or 1");out<<row[0]<<' ';}return out.str();});
        std::istringstream in(text);bits.assign(all_sites.size(),std::vector<int>(1));for(auto& row:bits)in>>row[0];
      }
      if(o.has("initial_statefile") || o.get("initial_state","zero")=="plus")for(std::size_t i=0;i<sites.size();++i) {
        auto shape=tcapi::shape(ctx,V[i]);std::vector<Elem> data(2);
        if(o.has("initial_statefile"))data[bits[sites[i]][0]]=1;
        else data[0]=data[1]=1/std::sqrt(2.0);
        V[i]=tcapi::assign_from_range<Tensor>(ctx,shape,data.begin(),[&shape](const auto& c){return c[shape.size()-1];});
      }
    }
    auto measurements=ite::measurement_ops(ctx,plan);auto actions=tnbp::MakeTrotterActions(plan,cfg.settings.trotter_order);
    std::vector<std::vector<std::vector<int>>> layer_batches;
    std::vector<std::vector<bool>> layer_touches;
    if(rank==0)ite::log(true,"compression: leaf-batches-v1");
    for(std::size_t l=0;l<plan.layers.size();++l) {
      layer_batches.push_back(tnbp::MakePauliLayerTruncationBatches(plan,l));
      layer_touches.emplace_back(all_sites.size(),false);
      for(auto e:plan.layers[l].edge_ids) {
        auto [a,b]=edges[e];layer_touches.back()[a]=layer_touches.back()[b]=true;
      }
      if(rank==0)for(std::size_t b=0;b<layer_batches.back().size();++b) {
        std::ostringstream ids;for(int e:layer_batches.back()[b])ids<<' '<<e;
        ite::log(true,"compression_batch layer=",l," batch=",b," edges=",ids.str());
      }
    }
    auto bp_order=tnbp::MakeBPCommunicationOrder(edges,owners,rank);
    ite::Window window{next};std::optional<std::size_t> last_summary;
    // A stage-local cache is reused for forward/reverse actions and all its steps.
    std::map<std::pair<std::size_t,double>,tnbp::ScaledPauliTPO<Tensor>> cache;double cached_dt=-1;
    for(;next<end;++next) {
      const double dt=ite::step_size(cfg.schedule,next);if(dt!=cached_dt){cache.clear();cached_dt=dt;}
      for(auto action:actions) {
        const auto key=std::make_pair(action.layer_id,action.time_fraction);
        auto it=cache.find(key);if(it==cache.end())it=cache.emplace(key,tnbp::MakePauliLayerExponentialTPO<Tensor>(ctx,plan,action.layer_id,dt*action.time_fraction,sites)).first;
        tnbp::AttachTPO(ctx,edges,it->second.tensors,V,sites,owners,E,edge_idx,comm);
        operator_log_scale+=it->second.log_scale;ite::require(std::isfinite(operator_log_scale),"operator log scale overflow");
        auto converge_bp=[&](std::size_t batch_index) {
          if(edges.empty()) return;
          bool converged=false;Real residual=std::numeric_limits<Real>::infinity();
          for(std::size_t k=0;k<cfg.settings.max_bp_iterations;++k) {
            tnbp::UpdateBPMessages(ctx,edges,V,sites,E,edge_idx,comm,F);
            tnbp::ExchangeBPMessages(ctx,edges,owners,edge_idx,bp_order,comm,F);
            // MPI sends tensor values, not backend-internal scale factors.
            // Rebuild both copies in the same value representation: untouched
            // edges no longer get canonical diagonal messages from Truncation.
            for(auto& message:F) {
              auto shape=tcapi::shape(ctx,message);
              std::vector<Elem> values(tcapi::size(ctx,message));
              tcapi::to_range(ctx,message,values.begin(),[&shape](const auto& c){return tnbp::address_from_coor(shape,c);});
              message=tcapi::assign_from_range<Tensor>(ctx,shape,values.begin(),[&shape](const auto& c){return tnbp::address_from_coor(shape,c);});
            }
            tnbp::BeliefPropagationDifference(ctx,edges,owners,E,F,edge_idx,comm,residual);E.swap(F);
            if(residual<cfg.settings.bp_tolerance){converged=true;break;}
          }
          if(!converged && rank==0)ite::log(false,"WARNING BP_NOT_CONVERGED step=",next+1," layer=",action.layer_id," before_batch=",batch_index," residual=",ite::number(residual));
        };
        converge_bp(0);
        const auto& batches=layer_batches[action.layer_id];
        for(std::size_t b=0;b<batches.size();++b) {
          if(b) converge_bp(b);
          auto result=tnbp::TruncationBatch(ctx,edges,batches[b],V,sites,owners,E,edge_idx,comm,
            static_cast<BondDim>(cfg.settings.max_bond_dim),Real(cfg.settings.sv_min),Real(cfg.settings.truncation_error));
          for(std::size_t j=0;j<result.edge_ids.size();++j)
            if(owners.at(edges[result.edge_ids[j]].first)==rank)window.add(result.truncation_errors[j]);
        }
        // Single-site terms can act on vertices not visited by any batch.
        for(std::size_t i=0;i<sites.size();++i)if(!layer_touches[action.layer_id][sites[i]]) {
          auto norm=tcapi::normalize(ctx,V[i]);ite::require(std::isfinite(norm) && norm>0,"zero/nonfinite state");
        }
      }
      operator_log_scale-=dt*cfg.h.identity_coefficient;ite::require(std::isfinite(operator_log_scale),"identity log scale overflow");
      if((next+1)%measure_every==0) {
        ite::summary(ctx,cfg,plan,measurements,V,sites,owners,E,edge_idx,window,next+1,comm);last_summary=next+1;
      }
    }
    if(!last_summary || *last_summary!=next)ite::summary(ctx,cfg,plan,measurements,V,sites,owners,E,edge_idx,window,next,comm);
    if(o.has("portable_save")) {
      auto saved=ite::save(o.get("portable_save"),{"",plan.plan_digest,cfg,next,operator_log_scale},comm,[&](const std::string& target){
        return tnbp::SaveTPSPortable(ctx,target,edges,V,sites,owners,E,edge_idx,comm);
      });
      if(rank==0)ite::log(false,"SAVED step=",next," digest=",saved.state_digest," directory=",o.get("portable_save"));
    }
    if(rank==0)ite::log(false,"DONE step=",next," tau=",ite::elapsed_time(cfg.schedule,next)," operator_log_scale=",operator_log_scale);
  } catch(const std::exception& e) {
    std::cerr<<tnbp::make_timestamp()<<" ERROR rank="<<rank<<" "<<e.what()<<std::endl;MPI_Abort(comm,1);
  }
  MPI_Finalize();return 0;
}
