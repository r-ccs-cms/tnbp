#include "../common/device.h"
#include "tnbp/ptns/function/mcmc_sampling.h"
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
using Scalar=std::complex<double>;
#ifdef TNBP_TEST_CUDA
using Tensor=tcapi::cuda::Tensor<Scalar>;
#else
using Tensor=gqten::tensor<Scalar>;
#endif
// Count actual collectives, not logical audit events. PMPI keeps the normal
// implementation while detecting accidental synchronization in the hot path.
std::size_t allreduce_calls=0,bcast_calls=0,parent_collectives=0,comm_creations=0,comm_frees=0;
extern "C" int MPI_Comm_create_group(MPI_Comm parent,MPI_Group group,int tag,MPI_Comm* result) {
    ++comm_creations;return PMPI_Comm_create_group(parent,group,tag,result);
}
extern "C" int MPI_Comm_free(MPI_Comm* comm) {++comm_frees;return PMPI_Comm_free(comm);}
MPI_Comm monitored_parent=MPI_COMM_NULL;
extern "C" int MPI_Bcast(void* buffer,int count,MPI_Datatype type,int root,MPI_Comm comm) {
    ++bcast_calls;if(monitored_parent!=MPI_COMM_NULL && comm==monitored_parent)++parent_collectives;return PMPI_Bcast(buffer,count,type,root,comm);
}
extern "C" int MPI_Allreduce(const void* send,void* receive,int count,MPI_Datatype type,MPI_Op op,MPI_Comm comm) {
    ++allreduce_calls;if(monitored_parent!=MPI_COMM_NULL && comm==monitored_parent)++parent_collectives;return PMPI_Allreduce(send,receive,count,type,op,comm);
}
using Context=tcapi::context_handle_t<Tensor>;
using Edges=std::vector<std::pair<int,int>>;
using Labels=tcapi::List<tcapi::bond_label_t<Tensor>>;
#include "../common/dense_fixture.h"

// Enumerate virtual indices independently of all graph/label contractions.
Scalar exact(const Edges& edges,const std::vector<int>& state) {
    Scalar sum=0.;
    for(std::size_t bits=0;bits<(std::size_t(1)<<edges.size());++bits) {
        Scalar product=1.;
        for(std::size_t j=0;j<state.size();++j) {
            std::vector<int> bonds;
            for(std::size_t e=0;e<edges.size();++e)
                if(edges[e].first==int(j)||edges[e].second==int(j)) bonds.push_back((bits>>e)&1);
            product*=entry(j,state[j],bonds);
        }
        sum+=product;
    }
    return sum;
}
std::vector<Tensor> tensors(Context& ctx,const Edges& edges,const std::vector<int>& dims,
                             const std::vector<int>& ids) {
    std::vector<Tensor> result;
    for(int j:ids) {
        auto bonds=tnbp::GetSurroundingBondIndex(j,edges);
        tcapi::shape_t<Tensor> shape(bonds.size()+1,2);shape.back()=dims[j];
        auto t=tcapi::zeros<Tensor>(ctx,shape);
        for(int bits=0;bits<(1<<bonds.size());++bits) for(int s=0;s<dims[j];++s) {
            tcapi::elem_coors_t<Tensor> c;std::vector<int> b;
            for(std::size_t k=0;k<bonds.size();++k) { c.push_back((bits>>k)&1);b.push_back((bits>>k)&1); }
            c.push_back(s);tcapi::set_elem(ctx,t,c,entry(j,s,b));
        }
        result.push_back(std::move(t));
    }
    return result;
}
void require(bool ok,const char* message) {if(!ok) throw std::runtime_error(message);}

// Compile the actual commit function with a noncopyable value. Pointer identity
// proves that commit consumes the prepared slice rather than rebuilding it.
void commit_move_only_check(MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    const std::vector<int> line{2,0,3},selected{1,2,0},dimensions{2,3,1};
    for(bool empty_root:{false,true}) {
        std::map<int,int> owners;std::vector<int> ids,target_ids;
        for(int site=0;site<4;++site) {
            owners[site]=empty_root ? size-1:(site+1)%size;
            if(owners[site]==rank) {
                ids.push_back(site);
                if(site!=1) target_ids.push_back(site);
            }
        }
        std::reverse(ids.begin(),ids.end());
        std::vector<int> states(ids.size(),-7);
        std::vector<std::vector<std::unique_ptr<int>>> slices(line.size());
        std::vector<std::unique_ptr<int>> target(target_ids.size());
        std::vector<const int*> selected_storage(line.size(),nullptr);
        for(std::size_t j=0;j<line.size();++j) if(owners.at(line[j])==rank) {
            for(int s=0;s<dimensions[j];++s) slices[j].push_back(std::make_unique<int>(100*line[j]+s));
            selected_storage[j]=slices[j][selected[j]].get();
        }
        tnbp::mcmc_commit_window_state(slices,ids,owners,line,selected,rank,target,target_ids,states);
        for(std::size_t j=0;j<line.size();++j) if(owners.at(line[j])==rank) {
            const auto ti=std::find(target_ids.begin(),target_ids.end(),line[j])-target_ids.begin();
            const auto wi=std::find(ids.begin(),ids.end(),line[j])-ids.begin();
            require(target[ti].get()==selected_storage[j] && *target[ti]==100*line[j]+selected[j],"commit did not reuse selected raw slice storage");
            require(!slices[j][selected[j]] && states[wi]==selected[j],"commit did not move slice/update owner state");
            for(int s=0;s<dimensions[j];++s) if(s!=selected[j])
                require(slices[j][s] && *slices[j][s]==100*line[j]+s,"commit touched an unselected raw slice");
        }
        for(std::size_t wi=0;wi<ids.size();++wi) if(ids[wi]==1)
            require(states[wi]==-7,"commit modified an off-line state");
    }
    if(rank==0) std::cout<<"move-only commit PASS: original storage reused, selected slices consumed, owner/state mapping and unselected slices preserved\n";
}

// Enumerate all configurations and all windows on a genuinely distributed graph.
void conditional_check(Context& ctx,const Edges& edges,const std::vector<int>& dims,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    std::vector<int> ids,line{0,1,2,3};std::map<int,int> owners;
    for(auto j:line) {owners[j]=j%size;if(owners[j]==rank) ids.push_back(j);}
    auto W=tensors(ctx,edges,dims,ids);
    std::vector<std::vector<int>> groups{{0},{1},{2},{3}};
    const std::vector<Tensor> empty;const std::vector<int> ie;
    std::vector<std::vector<Tensor>> slices;std::vector<std::vector<int>> axes;
    std::vector<int> physical_dims,group_owners;std::vector<tcapi::bond_dim_t<Tensor>> bond_dims;
    tnbp::mcmc_prepare_heat_bath_cache(ctx,edges,W,ids,owners,empty,ie,empty,ie,
        line,groups,std::vector<int>(ids.size(),0),comm,slices,axes,physical_dims,bond_dims);
    for(auto j:line) {
        group_owners.push_back(owners.at(j));
        require(slices[j].size()==std::size_t(rank==owners.at(j) ? dims[j]:0),"group cache not owner-local");
    }
    int count=1;for(int d:dims) count*=d;
    double max_error=0;
    for(int configuration=0;configuration<count;++configuration) {
        auto x=configuration;std::vector<int> state(4);
        for(int j=0;j<4;++j) {state[j]=x%dims[j];x/=dims[j];}
        for(int direction:{0,1}) for(int width:{1,2}) {
            std::vector<std::size_t> order{0,1,2,3};if(direction) std::reverse(order.begin(),order.end());
            std::vector<std::vector<Tensor>> suffix;std::vector<std::vector<int>> suffix_axes;
            std::vector<std::size_t> audit(10,0);
            tnbp::mcmc_build_heat_bath_environments(ctx,slices,axes,group_owners,state,order,bond_dims,
                comm,suffix,suffix_axes,width,&audit);
            require(audit[0]==std::size_t(3-width) && audit[5]==1,"opposite environments not built exactly once");
            if(width==2) require(suffix[1].empty() && suffix_axes[1].empty(),"unused pair environment constructed");
            std::vector<Scalar> fingerprints;
            for(const auto& cache:suffix) for(const auto& t:cache) {
                auto shape=tcapi::shape(ctx,t);std::size_t n=1;for(auto d:shape) n*=d;
                for(std::size_t k=0;k<n;++k) {
                    auto index=k;tcapi::elem_coors_t<Tensor> coors;
                    for(auto d:shape) {coors.push_back(index%d);index/=d;}
                    fingerprints.push_back(tcapi::get_elem(ctx,t,coors));
                }
            }
            std::vector<Tensor> prefix;std::vector<int> prefix_axes;
            for(int pos=0;pos+width<=4;++pos) {
                auto i=order[pos],j=width==2 ? order[pos+1]:i;auto end=std::size_t(pos+width);
                std::vector<double> a,w;std::vector<std::vector<Tensor>> candidates;std::vector<int> candidate_axes;
                tnbp::mcmc_evaluate_window_weights(ctx,slices[i],axes[i],group_owners[i],dims[i],
                    width==2 ? slices[j]:empty,width==2 ? axes[j]:ie,group_owners[j],width==2 ? dims[j]:1,width==2,
                    prefix,prefix_axes,pos>0,pos>0 ? group_owners[order[pos-1]]:group_owners[i],
                    suffix[end],suffix_axes[end],end<4,end<4 ? group_owners[order[end]]:group_owners[i],
                    bond_dims,comm,a,w,candidates,candidate_axes,&audit);
                require(a.size()==std::size_t(dims[i]*(width==2 ? dims[j]:1)),"candidate count");
                require(audit[0]==std::size_t(3-width) && audit[5]==1,"candidate rebuilt an exterior environment");
                require(audit[3]<=2,"candidate join included more than two cached regions");
                for(int c=0;c<dims[i];++c)
                    require(candidates[c].size()==std::size_t(rank==owners.at(i)),"prefix cache not owner-local");
                std::vector<double> reference(a.size());auto candidate=state;
                double total=0,ref_total=0;
                for(std::size_t k=0;k<a.size();++k) {
                    candidate[i]=k%dims[i];if(width==2) candidate[j]=k/dims[i];
                    const auto expected=std::abs(exact(edges,candidate));
                    max_error=std::max(max_error,std::abs(a[k]-expected)/(1+expected));
                    reference[k]=expected*expected;ref_total+=reference[k];total+=w[k];
                }
                for(std::size_t k=0;k<a.size();++k)
                    require(std::abs(w[k]/total-reference[k]/ref_total)<1.e-11,"conditional probability mismatch");
                const auto* pointer=candidates[state[i]].data();
                tnbp::mcmc_absorb_group(static_cast<std::size_t>(state[i]),candidates,candidate_axes,prefix,prefix_axes,&audit);
                require(prefix.data()==pointer,"candidate storage was copied rather than reused");
            }
            std::vector<Scalar> after;
            for(const auto& cache:suffix) for(const auto& t:cache) {
                auto shape=tcapi::shape(ctx,t);std::size_t n=1;for(auto d:shape) n*=d;
                for(std::size_t k=0;k<n;++k) {
                    auto index=k;tcapi::elem_coors_t<Tensor> coors;
                    for(auto d:shape) {coors.push_back(index%d);index/=d;}
                    after.push_back(tcapi::get_elem(ctx,t,coors));
                }
            }
            require(after==fingerprints,"candidate evaluation mutated an opposite environment");
            require(audit[4]==std::size_t(5-width),"intermediate reuse count");
        }
    }
    require(max_error<1.e-11,"candidate amplitude mismatch");
    if(rank==0) std::cout<<"distributed all configurations/candidates, both directions, widths 1/2: max relative amplitude error "<<max_error<<" PASS\n";
}

// Independent exact heat-bath oracle. Same specified roulette random draw, but
// independently enumerated amplitudes and probabilities; checks every window
// including overlapping pairs over repeated forward/backward sweeps.
void oracle_sweep(const Edges& edges,const std::vector<int>& dims,std::vector<int>& state,
                  std::vector<std::mt19937>& rngs,const std::map<int,int>& owners,int direction,int width,const std::vector<int>& line={0,1,2,3}) {
    auto order=line;if(direction) std::reverse(order.begin(),order.end());
    for(int pos=0;pos+width<=static_cast<int>(order.size());++pos) {
        auto i=order[pos],j=width==2 ? order[pos+1]:i;auto candidate=state;
        std::vector<long double> w(dims[i]*(width==2 ? dims[j]:1));long double total=0;
        for(std::size_t k=0;k<w.size();++k) {
            candidate[i]=k%dims[i];if(width==2) candidate[j]=k/dims[i];
            w[k]=std::norm(exact(edges,candidate));total+=w[k];
        }
        auto draw=std::generate_canonical<long double,std::numeric_limits<long double>::digits>(rngs.at(owners.at(i)))*total;
        std::size_t selected=w.size()-1;long double cumulative=0;
        for(std::size_t k=0;k<w.size();++k) {cumulative+=w[k];if(draw<cumulative) {selected=k;break;}}
        state[i]=selected%dims[i];if(width==2) state[j]=selected/dims[i];
    }
}
void distributed_check(Context& ctx,const Edges& edges,const std::vector<int>& dims,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    for(bool empty_root:{false,true}) for(int width:{1,2}) for(bool middle:{false,true}) {
        std::vector<int> ids,line{0,1,2,3};std::map<int,int> owners;
        for(int j:line) {owners[j]=empty_root ? size-1 : j%size;if(owners[j]==rank) ids.push_back(j);}
        auto W=tensors(ctx,edges,dims,ids);std::vector<int> state(ids.size()),expected(4);
        std::vector<Tensor> target;for(auto& t:W) {
            Labels l(tcapi::order(ctx,t));std::iota(l.begin(),l.end(),0);Tensor attached;Labels la;
            tnbp::mcmc_attach_window_state(ctx,t,l,l.back(),0,attached,la);target.push_back(std::move(attached));
        }
        std::vector<Tensor> upper,lower;std::vector<int> upper_ids,lower_ids;
        if(middle) for(std::size_t k=0;k<ids.size();++k) {
            if(ids[k]==0) {lower.push_back(tcapi::copy(ctx,target[k]));lower_ids.push_back(0);}
            if(ids[k]==3) {upper.push_back(tcapi::copy(ctx,target[k]));upper_ids.push_back(3);}
        }
        const std::vector<int> update_line=middle ? std::vector<int>{1,2}:line;
        const std::vector<std::vector<int>> groups=middle ? std::vector<std::vector<int>>{{1,0},{2,3}}:
            std::vector<std::vector<int>>{{0},{1},{2},{3}};
        std::mt19937 rng(913);std::vector<std::mt19937> oracle_rngs(size,std::mt19937(913));
        for(int step=0;step<32;++step) {
            auto direction=step%2;
            oracle_sweep(edges,dims,expected,oracle_rngs,owners,direction,width,update_line);
            std::vector<std::size_t> audit;
            tnbp::mcmc_update_line_heat_bath(ctx,edges,W,ids,owners,upper,upper_ids,target,ids,lower,lower_ids,
                update_line,groups,state,comm,rng,direction,width,&audit);
            const bool participating=std::any_of(owners.begin(),owners.end(),[&](const auto& item){return item.second==rank;});
            if(participating) {
            require(audit[0]==(middle ? std::size_t(width==1):update_line.size()-width-1) && audit[2]==update_line.size() && audit[5]==1,
                "exterior/group construction repeated inside sweep");
            require(audit[4]<=update_line.size()-width+1 && audit[3]<=4,"window reuse/region count");
            } else require(std::all_of(audit.begin(),audit.end(),[](auto n){return n==0;}),"nonparticipant executed line work");
            for(std::size_t k=0;k<ids.size();++k) {
                require(state[k]==expected[ids[k]],"MPI sweep differs from independent oracle");
                auto shape=tcapi::shape(ctx,target[k]);std::size_t elements=1;for(auto d:shape) elements*=d;
                for(std::size_t v=0;v<elements;++v) {
                    auto bits=v;tcapi::elem_coors_t<Tensor> c;std::vector<int> b;
                    for(auto d:shape) {c.push_back(bits%d);b.push_back(bits%d);bits/=d;}
                    require(std::abs(tcapi::get_elem(ctx,target[k],c)-entry(ids[k],state[k],b))<1.e-12,"owner tensor/state mismatch");
                }
            }
        }
    }
    if(rank==0) std::cout<<"MPI owner commit/overlapping pairs/empty-root oracle PASS\n";
}
void roulette_check(Context& ctx) {
    std::mt19937 rng(901);std::vector<double> w;
    for(double scale:{1.e-250,1.,1.e250}) {
        tnbp::mcmc_heat_bath_weights(std::vector<double>{0.,scale,2*scale},w);
        require(w==std::vector<double>({0.,.25,1.}),"amplitude scaling");
    }
    for(auto bad:std::vector<std::vector<double>>{{},{0.,0.},{-1.,2.},{1.,std::numeric_limits<double>::infinity()}}) {
        bool caught=false;try {tnbp::mcmc_draw_heat_bath(bad,rng);}catch(const std::runtime_error&){caught=true;}
        require(caught,"invalid weight accepted");
    }
    int count=0;for(int k=0;k<20000;++k) {auto s=tnbp::mcmc_draw_heat_bath(std::vector<double>{0.,1.,3.},rng);require(s!=0,"zero bin selected");count+=s==2;}
    require(std::abs(count/20000.-.75)<.015,"roulette frequencies");
    std::cout<<"roulette/scaling/invalid weights PASS\n";
}


void support_check(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    const Edges edge{{0,1}};const std::vector<int> line{0,1};
    const std::vector<std::vector<int>> groups{{0},{1}};
    const std::vector<Tensor> empty;const std::vector<int> ie;
    // Both distinct owners and an entirely empty rank zero.
    for(bool empty_root:{false,true}) {
        const std::map<int,int> owners{{0,empty_root ? size-1:0},{1,size-1}};
        std::vector<int> ids;for(auto site:line) if(owners.at(site)==rank) ids.push_back(site);
        std::vector<Tensor> W;
        for(auto site:ids) {
            auto t=tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{2,2});
            for(int s=0;s<2;++s)
                tcapi::set_elem(ctx,t,tcapi::elem_coors_t<Tensor>{site==0 ? s:1-s,s},Scalar(1));
            W.push_back(std::move(t));
        }
        for(int width:{1,2}) {
            std::vector<int> states;std::vector<Tensor> target;
            for(auto site:ids) {states.push_back(site);target.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{2}));}
            std::mt19937 rng(713+rank);bool saw01=false,saw10=false;
            for(int step=0;step<64;++step) {
                tnbp::mcmc_update_line_heat_bath(ctx,edge,W,ids,owners,empty,ie,target,ids,empty,ie,
                    line,groups,states,comm,rng,step%2,width);
                int local[2]={0,0},selected[2]={0,0};
                for(std::size_t k=0;k<ids.size();++k) local[ids[k]]=states[k];
                MPI_Allreduce(local,selected,2,MPI_INT,MPI_SUM,comm);
                require(selected[0]+selected[1]==1,"pair left fixed-charge support");
                if(width==1) require(selected[0]==0 && selected[1]==1,"one-site fixed-charge state not frozen");
                saw01=saw01 || selected[0]==0;saw10=saw10 || selected[0]==1;
            }
            if(width==2) require(saw01 && saw10,"two-site failed to transfer charge");
        }
        // All d=16 x d=16 candidates, with different charges, stay eligible.
        W.clear();std::vector<int> state(ids.size(),0);
        for(auto site:ids) {
            auto t=tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{1,16});
            for(int s=0;s<16;++s) tcapi::set_elem(ctx,t,tcapi::elem_coors_t<Tensor>{0,s},Scalar(site==0 ? s+1:2*s+1));
            W.push_back(std::move(t));
        }
        std::vector<std::vector<Tensor>> slices;std::vector<std::vector<int>> axes;
        std::vector<int> dims;std::vector<tcapi::bond_dim_t<Tensor>> bond_dims;
        tnbp::mcmc_prepare_heat_bath_cache(ctx,edge,W,ids,owners,empty,ie,empty,ie,line,groups,state,
            comm,slices,axes,dims,bond_dims);
        std::vector<double> a,w;std::vector<std::vector<Tensor>> candidates;std::vector<int> candidate_axes;
        tnbp::mcmc_evaluate_window_weights(ctx,slices[0],axes[0],owners.at(0),16,
            slices[1],axes[1],owners.at(1),16,true,
            empty,ie,false,owners.at(0),empty,ie,false,owners.at(1),bond_dims,comm,a,w,candidates,candidate_axes);
        require(a.size()==256 && candidates.size()==16,"p16 candidate/cache count");
        for(int t=0;t<16;++t) for(int s=0;s<16;++s)
            require(std::abs(a[s+16*t]-(s+1)*(2*t+1))<1.e-12 && w[s+16*t]>0,"p16 candidate filtered");
        if(empty_root && size>1) require(rank!=0 || (slices[0].empty() && slices[1].empty() && candidates[0].empty()),"empty root acquired tensors");
    }
    // A singleton line is skipped in pair mode and updated normally in one-site mode.
    std::map<int,int> owner{{0,size-1}};std::vector<int> ids,state;
    std::vector<Tensor> W,target;
    if(rank==size-1) {
        ids={0};state={0};W.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{1}));
        tcapi::set_elem(ctx,W.back(),tcapi::elem_coors_t<Tensor>{0},Scalar(1));
        target.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{}));
    }
    std::mt19937 rng(91);auto before=rng;
    tnbp::mcmc_update_line_heat_bath(ctx,Edges{},W,ids,owner,empty,ie,target,ids,empty,ie,
        std::vector<int>{0},std::vector<std::vector<int>>{{0}},state,comm,rng,0,2);
    require(rng==before,"singleton pair mode consumed random draw");
    tnbp::mcmc_update_line_heat_bath(ctx,Edges{},W,ids,owner,empty,ie,target,ids,empty,ie,
        std::vector<int>{0},std::vector<std::vector<int>>{{0}},state,comm,rng,1,1);
    if(rank==size-1) require(state==std::vector<int>{0} &&
        tcapi::get_elem(ctx,target[0],tcapi::elem_coors_t<Tensor>{})==Scalar(1),"singleton one-site commit");
    if(rank==0) std::cout<<"distributed fixed-charge support and all 256 p16 candidates PASS\n";
}

// Operation counts must depend on line length only during environment setup,
// even when the number of one-/two-site candidates changes from 1 to 3 or 9.
void cache_reuse_check(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    constexpr int n=8;
    Edges edges;std::vector<int> line;std::vector<std::vector<int>> groups;
    std::map<int,int> owners;std::vector<int> ids;
    for(int j=0;j<n;++j) {
        if(j) edges.emplace_back(j-1,j);
        line.push_back(j);groups.push_back({j});owners[j]=(j+1)%size;
        if(owners[j]==rank) ids.push_back(j);
    }
    const std::vector<Tensor> empty;const std::vector<int> ie;
    for(int d:{1,3}) for(int width:{1,2}) for(int direction:{0,1}) {
        std::vector<Tensor> W,target;std::vector<int> states(ids.size(),0);
        for(int j:ids) {
            auto shape=tcapi::shape_t<Tensor>(tnbp::GetSurroundingBondIndex(j,edges).size()+1,1);
            shape.back()=d;auto t=tcapi::zeros<Tensor>(ctx,shape);
            for(int s=0;s<d;++s) {
                tcapi::elem_coors_t<Tensor> c(shape.size(),0);c.back()=s;
                tcapi::set_elem(ctx,t,c,Scalar(1.+.1*j+.2*s));
            }
            W.push_back(std::move(t));shape.pop_back();target.push_back(tcapi::zeros<Tensor>(ctx,shape));
        }
        std::mt19937 rng(61+rank);std::vector<std::size_t> audit;
        const auto before_collectives=allreduce_calls,before_bcasts=bcast_calls;
        monitored_parent=comm;parent_collectives=0;
        tnbp::mcmc_update_line_heat_bath(ctx,edges,W,ids,owners,empty,ie,target,ids,empty,ie,
            line,groups,states,comm,rng,direction,width,&audit);
        require(parent_collectives==0,"line update used parent communicator collective");monitored_parent=MPI_COMM_NULL;
        const auto plan=tnbp::mcmc_build_edge_line_plan(edges,line,groups,owners,direction,width);
        std::size_t env_calls=0,window_calls=0,windows=0,metadata=0,joins=0;
        auto member=[&](const auto& ranks){return std::binary_search(ranks.begin(),ranks.end(),rank);};
        for(const auto& env:std::get<1>(plan))if(member(std::get<3>(env))) {
            env_calls+=std::get<2>(env).size();metadata+=std::get<2>(env).size();
        }
        for(const auto& win:std::get<2>(plan))if(member(std::get<3>(win))) {
            ++windows;
            for(const auto* region:{&std::get<0>(win),&std::get<1>(win)}) {
                const auto common=std::get<6>(*region),total=std::get<2>(*region).size();
                window_calls+=common+(total-common)*d;metadata+=total;
            }
            if(!std::get<2>(std::get<2>(win)).empty()){++joins;window_calls+=width==1?d:d*d;}
        }
        require(allreduce_calls==before_collectives,"native line added a shape/validation Allreduce");
        require(bcast_calls-before_bcasts==windows*std::size_t(width+2),"unexpected window collective");
        require(audit[0]==env_calls && audit[5]==1 && audit[2]==n,"exterior construction repeated");
        require(audit[1]==window_calls,"native window contractions differ from compiled plan");
        require(audit[3]<=2 && audit[4]==windows,"native region inputs/selected intermediate reuse");
        require(audit[6]==0 && audit[9]==env_calls+window_calls,"identity or extra native contraction");
        require(audit[8]==0,"single-node groups unexpectedly recomputed common subtrees");
        require(audit[13]==joins && audit[12]==metadata,"metadata rebuilt per candidate");
        std::size_t consumed=0;const auto& order=std::get<0>(plan);
        for(std::size_t pos=0;pos+width<=line.size();++pos) {
            if(pos && owners.at(line[order[pos-1]])==rank)++consumed;
            const auto end=pos+width;
            if(width==2) {
                if(end<line.size() && owners.at(line[order[end]])==rank)++consumed;
                if(owners.at(line[order[pos+1]])==rank)consumed+=d;
            }else if(end<line.size() && owners.at(line[order[end]])==rank)++consumed;
        }
        require(audit[11]==consumed,"last-use prefix/suffix/RHS was copied instead of moved");
        for(std::size_t k=0;k<ids.size();++k) {
            require(states[k]>=0 && states[k]<d,"cached product state bounds");
            require(std::abs(tcapi::get_elem(ctx,target[k],tcapi::elem_coors_t<Tensor>(tcapi::order(ctx,target[k]),0))-
                Scalar(1.+.1*ids[k]+.2*states[k]))<1.e-12,"cached product tensor/state mismatch");
        }
    }
    if(rank==0) std::cout<<"8-group native audit PASS: compiled environment/window counts, zero parent collectives and line Allreduces, selected intermediates reused, both directions\n";
}

// Noncontiguous parent owners {0,2}; rank 0 owns only the boundary, rank 2
// owns targets, and rank 1 waits outside the update until it has completed.
void line_communicator_check(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);if(size!=3)return;
    const Edges edges{{0,1},{1,2},{2,3}};
    const std::map<int,int> owners{{0,0},{1,2},{2,2},{3,0},{4,1}};
    const std::vector<int> dims{1,2,2,1,2},line{1,2};
    const std::vector<std::vector<int>> groups{{1,0},{2,3}};
    const auto created_before=comm_creations,freed_before=comm_frees;
    std::map<std::vector<int>,MPI_Comm> handles;
    tnbp::mcmc_acquire_region_communicators(comm,0,{{0,2}},handles);
    const auto created=comm_creations;
    tnbp::mcmc_acquire_region_communicators(comm,0,{{0,2}},handles);
    require(comm_creations==created,"duplicate registry creation");
    for(int width:{1,2})for(int direction:{0,1}) {
        std::vector<int> ids,upper_ids,lower_ids;
        for(const auto& item:owners)if(item.second==rank)ids.push_back(item.first);
        auto W=tensors(ctx,edges,dims,ids);std::vector<Tensor> target,upper,lower;
        std::vector<int> states(ids.size(),0),expected{0,0,0,0,1};
        for(std::size_t k=0;k<ids.size();++k) {
            if(ids[k]==4)states[k]=1;
            Labels labels(tcapi::order(ctx,W[k]));std::iota(labels.begin(),labels.end(),0);Labels result_labels;Tensor projected;
            tnbp::mcmc_attach_window_state(ctx,W[k],labels,labels.back(),states[k],projected,result_labels);
            if(ids[k]==0){upper_ids.push_back(0);upper.push_back(tcapi::copy(ctx,projected));}
            if(ids[k]==3){lower_ids.push_back(3);lower.push_back(tcapi::copy(ctx,projected));}
            target.push_back(std::move(projected));
        }
        std::mt19937 rng(901+rank);const auto rng_before=rng;
        std::vector<std::mt19937> oracle;for(int r=0;r<size;++r)oracle.emplace_back(901+r);
        oracle_sweep(edges,dims,expected,oracle,owners,direction,width,line);
        const auto original_edges=edges;const auto original_owners=owners;
        const auto reductions=allreduce_calls,broadcasts=bcast_calls;parent_collectives=0;monitored_parent=comm;
        int ready=1;if(rank==1)MPI_Recv(&ready,1,MPI_INT,2,92,comm,MPI_STATUS_IGNORE);
        std::vector<std::size_t> audit;
        tnbp::mcmc_update_line_heat_bath(ctx,edges,W,ids,owners,upper,upper_ids,target,ids,lower,lower_ids,
            line,groups,states,comm,rng,direction,width,&audit,&handles);
        if(rank==2)MPI_Send(&ready,1,MPI_INT,1,92,comm);
        require(parent_collectives==0,"parent collective inside scoped line update");monitored_parent=MPI_COMM_NULL;
        require(edges==original_edges && owners==original_owners,"replicated graph/owners changed");
        if(rank==1) {
            require(rng==rng_before && states==std::vector<int>{1},"nonparticipant state/RNG changed");
            require(allreduce_calls==reductions && bcast_calls==broadcasts && audit[9]==0,"nonparticipant communicated");
            require(tcapi::get_elem(ctx,target[0],tcapi::elem_coors_t<Tensor>{})==entry(4,1,{}),"unrelated tensor changed");
        } else {
            require(allreduce_calls==reductions && audit[5]==1,"region shape sharing or environments changed");
            for(std::size_t k=0;k<ids.size();++k)require(states[k]==expected[ids[k]],"remapped owner roulette differs from oracle");
            if(rank==0)require(rng==rng_before,"boundary-only owner drew RNG");
        }
    }
    require(comm_creations==created && comm_frees==freed_before,"line recreated/freed shared communicator");
    tnbp::mcmc_release_region_communicators(handles);
    require(comm_creations-created_before==std::size_t(rank!=1) && comm_frees-freed_before==std::size_t(rank!=1),"registry lifetime count");
    if(rank==0)std::cout<<"three-line communicator PASS: parent owners 0/2, boundary-only rank, rank1 delayed, unchanged global metadata and unrelated state, both widths/directions\n";
}

void review_fixes_check(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    for(bool empty_root:{false,true}) {
        const int oa=empty_root ? size-1:0,ob=size-1;
        // Closing one shared leg must allocate no identity, copy two inputs
        // once each, and execute only one graph contraction, independent of chi.
        for(int chi:{1,7,64}) {
            std::vector<Tensor> a,b,out;
            if(rank==oa) a.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{chi}));
            if(rank==ob) b.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{chi}));
            Scalar expected=0;
            for(int k=0;k<chi;++k) {
                auto av=entry(0,0,std::vector<int>{k}),bv=entry(1,1,std::vector<int>{k});
                expected+=av*bv;
                if(rank==oa) tcapi::set_elem(ctx,a[0],tcapi::elem_coors_t<Tensor>{k},av);
                if(rank==ob) tcapi::set_elem(ctx,b[0],tcapi::elem_coors_t<Tensor>{k},bv);
            }
            // Immutability compares stored values before/after, rather than a
            // separately optimized floating-point reevaluation of entry().
            std::vector<Scalar> before_a,before_b;
            double fixture_rounding=0.;
            for(int k=0;k<chi;++k) {
                if(rank==oa) {
                    before_a.push_back(tcapi::get_elem(ctx,a[0],tcapi::elem_coors_t<Tensor>{k}));
                    fixture_rounding=std::max(fixture_rounding,std::abs(before_a.back()-entry(0,0,std::vector<int>{k})));
                }
                if(rank==ob) {
                    before_b.push_back(tcapi::get_elem(ctx,b[0],tcapi::elem_coors_t<Tensor>{k}));
                    fixture_rounding=std::max(fixture_rounding,std::abs(before_b.back()-entry(1,1,std::vector<int>{k})));
                }
            }
            require(fixture_rounding<1.e-12,"stored fixture differs from reference");
            if(fixture_rounding>0.) std::cout<<"fixture reevaluation rounding rank "<<rank<<" chi "<<chi<<": "<<fixture_rounding<<" before contraction\n";
            std::vector<int> axes;std::vector<std::size_t> audit(10,0);
            tnbp::mcmc_contract_cached_regions(ctx,a,std::vector<int>{0},oa,b,std::vector<int>{0},ob,
                std::vector<tcapi::bond_dim_t<Tensor>>{chi},comm,out,axes,1,&audit);
            require(axes.empty() && audit[6]==0 && audit[7]==2 && audit[9]==1,"single-leg direct-edge work counts");
            if(rank==oa) require(std::abs(tcapi::get_elem(ctx,out[0],tcapi::elem_coors_t<Tensor>{})-expected)<1.e-11*(1+std::abs(expected)),"single-leg result");
            for(int k=0;k<chi;++k) {
                if(rank==oa) require(tcapi::get_elem(ctx,a[0],tcapi::elem_coors_t<Tensor>{k})==before_a[k],"left cache mutated");
                if(rank==ob) require(tcapi::get_elem(ctx,b[0],tcapi::elem_coors_t<Tensor>{k})==before_b[k],"right cache mutated");
            }
        }
        // Shared leg is last on both operands, after different external legs.
        std::vector<Tensor> a,b,out;
        if(rank==oa) {
            a.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{2,4}));
            for(int i=0;i<2;++i) for(int k=0;k<4;++k)
                tcapi::set_elem(ctx,a[0],tcapi::elem_coors_t<Tensor>{i,k},entry(0,i,std::vector<int>{k}));
        }
        if(rank==ob) {
            b.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{3,4}));
            for(int j=0;j<3;++j) for(int k=0;k<4;++k)
                tcapi::set_elem(ctx,b[0],tcapi::elem_coors_t<Tensor>{j,k},entry(1,j,std::vector<int>{k}));
        }
        std::vector<int> axes;std::vector<std::size_t> audit(10,0);
        tnbp::mcmc_contract_cached_regions(ctx,a,std::vector<int>{0,2},oa,b,std::vector<int>{1,2},ob,
            std::vector<tcapi::bond_dim_t<Tensor>>{2,3,4},comm,out,axes,1,&audit);
        require(axes==std::vector<int>({0,1}) && audit[6]==0 && audit[9]==1,"direct-edge external-axis order");
        if(rank==oa) for(int i=0;i<2;++i) for(int j=0;j<3;++j) {
            Scalar expected=0;for(int k=0;k<4;++k) expected+=entry(0,i,std::vector<int>{k})*entry(1,j,std::vector<int>{k});
            require(std::abs(tcapi::get_elem(ctx,out[0],tcapi::elem_coors_t<Tensor>{i,j})-expected)<1.e-11,"direct-edge components");
        }
    }
    // Physical leaves have labels/edges only, no owner or numerical tensor.
    const Edges edge{{0,1}};const std::vector<int> line{0,1},ie;
    const std::vector<std::vector<int>> groups{{0},{1}};const std::vector<Tensor> empty;
    std::map<int,int> owners{{0,0},{1,size-1}};std::vector<int> ids;
    for(int site:line) if(owners[site]==rank) ids.push_back(site);
    auto W=tensors(ctx,edge,std::vector<int>{2,3},ids);
    Edges graph;std::vector<Tensor> raw;std::vector<int> raw_ids,leaves,dims;std::map<int,int> ranks;
    tnbp::mcmc_prepare_heat_bath_graph(ctx,edge,W,ids,owners,empty,ie,empty,ie,line,groups,
        std::vector<int>(ids.size(),0),comm,graph,raw,raw_ids,ranks,leaves,dims);
    require(raw_ids==ids && raw.size()==ids.size() && ranks==owners,"unused numerical physical leaf still prepared");
    for(int leaf:leaves) require(!ranks.count(leaf) && std::find(raw_ids.begin(),raw_ids.end(),leaf)==raw_ids.end(),"physical leaf has data or owner");
    if(rank==0) std::cout<<"review fixes PASS: chi1/7/64 direct edges, axis order, immutable inputs and metadata-only leaves\n";
}

// Interleaved shared legs of unequal sizes expose reversed fusion order.
void planner_fusion_check(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    const std::vector<int> aa{0,1,3,5},bb{1,2,3,4,5};
    const std::vector<tcapi::bond_dim_t<Tensor>> dims{2,3,2,4,2,2};
    const int oa=size-1,ob=0;std::vector<Tensor> a,b,out;std::vector<int> axes;
    auto make=[&](const std::vector<int>& tokens,int id) {
        tcapi::shape_t<Tensor> shape;std::size_t count=1;
        for(auto t:tokens){shape.push_back(dims[t]);count*=dims[t];}
        auto tensor=tcapi::zeros<Tensor>(ctx,shape);
        for(std::size_t v=0;v<count;++v) {
            auto index=v;std::vector<int> indices;tcapi::elem_coors_t<Tensor> coors;
            for(auto d:shape){coors.push_back(index%d);indices.push_back(index%d);index/=d;}
            tcapi::set_elem(ctx,tensor,coors,entry(id,0,indices));
        }
        return tensor;
    };
    if(rank==oa)a.push_back(make(aa,0));if(rank==ob)b.push_back(make(bb,1));
    std::vector<std::vector<std::size_t>> trace;
    tnbp::mcmc_contract_cached_regions(ctx,a,aa,oa,b,bb,ob,dims,comm,out,axes,1,nullptr,&trace);
    require(axes==std::vector<int>({0,2,4}) && trace[0][2]==3 && trace[0][3]==192 && trace[0][4]==8,
        "prepared fused join metadata/cost mismatch");
    if(rank==oa) for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<2;++k) {
        Scalar expected=0.;
        for(int x=0;x<3;++x)for(int y=0;y<4;++y)for(int z=0;z<2;++z)
            expected+=entry(0,0,{i,x,y,z})*entry(1,0,{x,j,y,k,z});
        auto actual=tcapi::get_elem(ctx,out[0],tcapi::elem_coors_t<Tensor>{i,j,k});
        require(std::abs(actual-expected)<1.e-11*(1+std::abs(expected)),"shared-leg fusion changed index order");
    }
    if(rank==0)std::cout<<"prepared fusion PASS: shared dimensions 3/4/2, interleaved external legs, direct complex sum\n";
}

// Each fatal case runs in its own bounded MPI job; success is a test failure.
void fatal_failure_check(Context& ctx,MPI_Comm comm,int failure) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    const int owner=size-1;std::vector<Tensor> empty,target;
    std::vector<std::vector<Tensor>> output;std::vector<int> axes;
    if(failure>=5) {
        const Edges edges{{0,1}};const std::vector<int> line{0,1},none;
        const std::vector<std::vector<int>> groups{{0},{1}};
        std::map<int,int> owners{{0,0},{1,owner}};
        std::vector<int> ids,states;std::vector<Tensor> W,commit;
        for(int site:line) if(owners.at(site)==rank) {
            ids.push_back(site);states.push_back(failure==5 && site==1 ? 2:0);
            W.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{failure==6 && site==1 ? 3:2,2}));
            commit.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{2}));
        }
        if(failure==7 && rank==owner && !W.empty()) W.pop_back();
        if(failure==8 && rank==owner && !commit.empty()) commit.pop_back();
        std::mt19937 rng(1);
        tnbp::mcmc_update_line_heat_bath(ctx,edges,W,ids,owners,empty,none,commit,ids,empty,none,
            line,groups,states,comm,rng,0,1);
    } else if(failure==0) {
        if(rank==0)target.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{1}));
        std::vector<double> amplitudes,weights;
        tnbp::mcmc_evaluate_window_weights(ctx,target,std::vector<int>{0},0,1,empty,std::vector<int>{},0,1,false,
            empty,std::vector<int>{},false,0,empty,std::vector<int>{0},true,owner,
            std::vector<tcapi::bond_dim_t<Tensor>>{1},comm,amplitudes,weights,output,axes);
    } else if(failure==4) {
        if(rank==owner) for(int dim:{2,3}) target.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{dim}));
        tnbp::mcmc_prepare_group_candidates(ctx,target,std::vector<int>{0},owner,2,
            std::vector<std::vector<Tensor>>{},std::vector<std::vector<int>>{},std::vector<int>{},
            empty,std::vector<int>{},false,owner,std::vector<tcapi::bond_dim_t<Tensor>>{2},comm,output,axes);
    } else {
        std::vector<std::vector<Tensor>> boundary(2);
        if(rank==owner) target.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{failure==1?3:2}));
        if(rank==0) {
            boundary[0].push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{2,2}));
            boundary[1].push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{2}));
        }
        std::vector<std::size_t> audit(failure==2 && rank==0 ? 0:15,0);
        std::vector<tcapi::bond_dim_t<Tensor>> dims{2,2};if(failure==3 && rank==owner)dims[0]=0;
        tnbp::mcmc_prepare_group_candidates(ctx,target,std::vector<int>{0},owner,1,
            boundary,std::vector<std::vector<int>>{{0,1},{1}},std::vector<int>{0,0},
            empty,std::vector<int>{},false,owner,dims,comm,output,axes,1,&audit);
    }
}

void scalar_failure_check(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);const int owner=size-1;
    std::vector<Tensor> first,empty;
    if(rank==owner) for(double value:{1.,std::numeric_limits<double>::quiet_NaN(),2.}) {
        auto t=tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{});
        tcapi::set_elem(ctx,t,tcapi::elem_coors_t<Tensor>{},Scalar(value));first.push_back(std::move(t));
    }
    bool caught=false;const auto before=allreduce_calls,bcasts_before=bcast_calls;
    std::vector<double> amplitudes{17},weights{19};std::vector<std::vector<Tensor>> candidates;std::vector<int> axes;
    try {
        tnbp::mcmc_evaluate_window_weights(ctx,first,std::vector<int>{},owner,3,empty,std::vector<int>{},0,1,false,
            empty,std::vector<int>{},false,0,empty,std::vector<int>{},false,0,
            std::vector<tcapi::bond_dim_t<Tensor>>{},comm,amplitudes,weights,candidates,axes);
    } catch(const std::runtime_error&) {caught=true;}
    require(caught && allreduce_calls==before && bcast_calls-bcasts_before==1 && amplitudes==std::vector<double>{17} && weights==std::vector<double>{19},
        "scalar failure should notify only at the window end without Allreduce");
    if(rank==0)std::cout<<"scalar failure PASS: terminal notification, unchanged outputs, no Allreduce\n";
}

// A three-column square strip, with D != chi. Reference amplitudes enumerate
// original edge indices directly; no contraction/cache helper is used.
void square_order_check(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    constexpr int D=2,chi=3;
    Edges edges{{0,1},{1,2},{3,4},{4,5},{6,7},{7,8},
                {0,3},{1,4},{2,5},{0,6},{1,7},{2,8}};
    std::vector<int> dimensions{D,D,chi,chi,chi,chi,D,D,D,D,D,D};
    std::vector<double> exact_amplitude(8,0.);
    std::size_t combinations=1;for(int d:dimensions) combinations*=d;
    std::vector<Scalar> exact_complex(8,Scalar(0));
    if(rank==0) for(std::size_t config=0;config<combinations;++config) {
        auto index=config;std::vector<std::vector<int>> bonds(9);
        for(std::size_t e=0;e<edges.size();++e) {
            const int value=index%dimensions[e];index/=dimensions[e];
            bonds[edges[e].first].push_back(value);bonds[edges[e].second].push_back(value);
        }
        Scalar outer=1.;for(int j=3;j<9;++j) outer*=entry(j,0,bonds[j]);
        for(int state=0;state<8;++state) {
            auto product=outer;for(int j=0;j<3;++j) product*=entry(j,(state>>j)&1,bonds[j]);
            exact_complex[state]+=product;
        }
    }
    if(rank==0) for(int state=0;state<8;++state) exact_amplitude[state]=std::abs(exact_complex[state]);
    MPI_Bcast(exact_amplitude.data(),8,MPI_DOUBLE,0,comm);
    for(bool empty_root:{false,true}) for(bool lower_first:{false,true}) {
        std::map<int,int> owners;
        for(int j=0;j<9;++j) owners[j]=empty_root ? size-1:(j+1+(size>2 ? j/3:0))%size;
        std::vector<int> ids,upper_ids,lower_ids;
        std::vector<Tensor> W,upper,lower,target;
        for(int j=0;j<9;++j) if(owners[j]==rank) {
            auto bonds=tnbp::GetSurroundingBondIndex(j,edges);
            tcapi::shape_t<Tensor> shape;std::size_t elements=1;
            for(auto e:bonds) {shape.push_back(dimensions[e]);elements*=dimensions[e];}
            if(j<3) shape.push_back(2);
            auto tensor=tcapi::zeros<Tensor>(ctx,shape);
            for(std::size_t v=0;v<elements;++v) for(int state=0;state<(j<3 ? 2:1);++state) {
                auto bits=v;std::vector<int> b;tcapi::elem_coors_t<Tensor> coors;
                for(auto e:bonds) {b.push_back(bits%dimensions[e]);coors.push_back(b.back());bits/=dimensions[e];}
                if(j<3) coors.push_back(state);
                tcapi::set_elem(ctx,tensor,coors,entry(j,state,b));
            }
            if(j<3) {
                ids.push_back(j);W.push_back(std::move(tensor));shape.pop_back();
                target.push_back(tcapi::zeros<Tensor>(ctx,shape));
            } else if(j<6) {upper_ids.push_back(j);upper.push_back(std::move(tensor));}
            else {lower_ids.push_back(j);lower.push_back(std::move(tensor));}
        }
        const std::vector<int> line{0,1,2};
        std::vector<std::vector<int>> groups;
        for(int j:line) groups.push_back(lower_first ? std::vector<int>{j,j+6,j+3}:std::vector<int>{j,j+3,j+6});
        std::vector<int> state(ids.size(),0),dims;
        std::vector<std::vector<Tensor>> slices;std::vector<std::vector<int>> axes;
        std::vector<std::vector<std::vector<Tensor>>> boundary;
        std::vector<std::vector<std::vector<int>>> boundary_axes;
        std::vector<std::vector<int>> boundary_owners;
        std::vector<tcapi::bond_dim_t<Tensor>> bond_dims;
        tnbp::mcmc_prepare_heat_bath_cache(ctx,edges,W,ids,owners,upper,upper_ids,lower,lower_ids,
            line,groups,state,comm,slices,axes,dims,bond_dims,nullptr,&boundary,&boundary_axes,&boundary_owners);
        // Exercise the actual candidate helper with a proper exterior side.
        // The middle target must execute boundary -> target -> boundary, in both
        // directions; pre-target work occurs once, not once per physical state.
        for(int direction:{0,1}) {
            std::vector<std::size_t> order{0,1,2};if(direction) std::reverse(order.begin(),order.end());
            std::vector<std::vector<Tensor>> suffix;std::vector<std::vector<int>> suffix_axes;
            const std::vector<int> selected(3,0),group_owners{owners[0],owners[1],owners[2]};
            tnbp::mcmc_build_heat_bath_environments(ctx,slices,axes,group_owners,selected,order,bond_dims,
                comm,suffix,suffix_axes,1,nullptr,&boundary,&boundary_axes,&boundary_owners);
            const auto side_index=order[2];
            std::vector<Tensor> three;
            if(rank==owners[1]) {
                for(int k:{0,1,0}) three.push_back(tcapi::copy(ctx,slices[1][k]));
            }
            std::vector<std::vector<std::uint64_t>> keys;
            std::vector<tnbp::mcmc_region_plan_arrays<int>> plans;
            std::vector<std::vector<Scalar>> reference_candidates;
            for(int count:{2,2,3}) {
            const bool cache_hit=count==2 && !keys.empty();
            std::vector<std::vector<Tensor>> candidates;std::vector<int> candidate_axes,steps;
            std::vector<std::vector<std::size_t>> trace;std::vector<std::size_t> audit(15,0);
            tnbp::mcmc_prepare_group_candidates(ctx,three,axes[1],owners[1],count,
                boundary[1],boundary_axes[1],boundary_owners[1],suffix[2],suffix_axes[2],true,owners[side_index],
                bond_dims,comm,candidates,candidate_axes,1,&audit,&steps,0,&trace,&keys,&plans);
            require(audit[10]==std::size_t(!cache_hit) && audit[14]==std::size_t(cache_hit) && audit[12]==3,
                "plan cache search/descriptor or per-step metadata reuse failed");
            if(count==2) {
            require(steps==std::vector<int>({0,-1,-1,-1,-1}),"square order is not common/target stage plan");
            require(trace.size()==5 && audit[6]==0 && audit[9]==5,"square order used extra contractions or identities");
            const std::size_t first_work=chi*chi*chi*D*D,target_work=chi*chi*D*D*D*D;
            for(std::size_t k=0;k<trace.size();++k) {
                require(trace[k][2]==std::size_t(k==0 ? 1:2),"square shared legs were not contracted together");
                require(trace[k][3]==(k==1 || k==3 ? target_work:first_work),"square MAC count differs from predicted powers");
                require(trace[k][4]==std::size_t(k==2 || k==4 ? chi*chi*D:chi*chi*D*D),"square intermediate exceeds expected size");
                const auto expected_owner=k==0 ? owners[side_index]:owners[1];
                require(trace[k][0]==std::size_t(expected_owner),"square intermediate owner convention");
            }
            } else {
                require(steps==std::vector<int>({0,0,-1,-1,-1}) && trace.size()==5,
                    "common subtree was not shared by all three states");
                std::size_t macs=0,peak=0;
                for(const auto& row:trace) {macs+=row[3];peak=std::max(peak,row[4]);}
                require(macs==756 && peak==72 && audit[7]==10,"weighted plan execution cost/copies mismatch");
                for(std::size_t k=0;k<trace.size();++k)
                    require(trace[k][0]==std::size_t(k<2 ? owners[side_index]:owners[1]),"common subtree owner mismatch");
            }
            for(const auto& candidate:candidates) require(candidate.size()==std::size_t(rank==owners[1]),"square candidate not target-local");
            if(rank==owners[1]) for(int k=0;k<count;++k) {
                const auto& t=candidates[k][0];const auto sh=tcapi::shape(ctx,t);
                std::size_t elements=1;for(auto d:sh) elements*=d;
                std::vector<Scalar> values;
                for(std::size_t v=0;v<elements;++v) {
                    auto index=v;tcapi::elem_coors_t<Tensor> coors;
                    for(auto d:sh) {coors.push_back(index%d);index/=d;}
                    values.push_back(tcapi::get_elem(ctx,t,coors));
                }
                if(count==2 && !cache_hit) reference_candidates.push_back(values);
                else for(std::size_t v=0;v<elements;++v)
                    require(std::abs(values[v]-reference_candidates[k%2][v])<1.e-11*(1+std::abs(values[v])),
                        "common-subtree plan changed candidate values");
            }
            }
        }
        // Compare every candidate in every window, in both directions, against
        // original-network sums for all eight initial configurations.
        for(int initial=0;initial<8;++initial) for(int direction:{0,1}) for(int width:{1,2}) {
            std::vector<int> selected(3);for(int j:line) selected[j]=(initial>>j)&1;
            const std::vector<int> group_owners{owners[0],owners[1],owners[2]};
            std::vector<std::size_t> order{0,1,2};if(direction) std::reverse(order.begin(),order.end());
            std::vector<std::vector<Tensor>> suffix;std::vector<std::vector<int>> suffix_axes;
            tnbp::mcmc_build_heat_bath_environments(ctx,slices,axes,group_owners,selected,order,bond_dims,
                comm,suffix,suffix_axes,width,nullptr,&boundary,&boundary_axes,&boundary_owners);
            std::vector<Tensor> prefix;std::vector<int> prefix_axes;
            for(std::size_t pos=0;pos+width<=3;++pos) {
                const auto a=order[pos],b=width==2 ? order[pos+1]:a,end=pos+width;
                std::vector<double> amplitudes,weights;std::vector<std::vector<Tensor>> candidates;std::vector<int> candidate_axes;
                tnbp::mcmc_evaluate_window_weights(ctx,slices[a],axes[a],owners[a],2,
                    slices[b],axes[b],owners[b],width==2 ? 2:1,width==2,
                    prefix,prefix_axes,pos>0,pos>0 ? owners[order[pos-1]]:owners[a],
                    suffix[end],suffix_axes[end],end<3,end<3 ? owners[order[end]]:owners[a],
                    bond_dims,comm,amplitudes,weights,candidates,candidate_axes,nullptr,
                    boundary[a],boundary_axes[a],boundary_owners[a],boundary[b],boundary_axes[b],boundary_owners[b]);
                for(std::size_t k=0;k<amplitudes.size();++k) {
                    auto candidate=selected;candidate[a]=k%2;if(width==2) candidate[b]=k/2;
                    const auto reference=exact_amplitude[candidate[0]+2*candidate[1]+4*candidate[2]];
                    require(std::abs(amplitudes[k]-reference)<1.e-11*(1+reference),"square full-state candidate amplitude differs from independent sum");
                }
                const auto* storage=candidates[selected[a]].data();
                tnbp::mcmc_absorb_group(std::size_t(selected[a]),candidates,candidate_axes,prefix,prefix_axes);
                require(prefix.data()==storage,"square selected intermediate not moved");
            }
        }
        // Repeated production sweeps, including overlapping pairs, must follow
        // an independent roulette oracle with one RNG stream per target owner.
        for(int width:{1,2}) {
            std::fill(state.begin(),state.end(),0);std::vector<int> expected(3,0);
            std::mt19937 rng(419+rank);std::vector<std::mt19937> oracle;
            for(int r=0;r<size;++r) oracle.emplace_back(419+r);
            for(int step=0;step<12;++step) {
                auto order=line;if(step%2) std::reverse(order.begin(),order.end());
                for(int pos=0;pos+width<=3;++pos) {
                    const int a=order[pos],b=width==2 ? order[pos+1]:a;
                    std::vector<double> w(width==2 ? 4:2);double total=0;
                    for(std::size_t k=0;k<w.size();++k) {
                        auto c=expected;c[a]=k%2;if(width==2) c[b]=k/2;
                        const auto amp=exact_amplitude[c[0]+2*c[1]+4*c[2]];w[k]=amp*amp;total+=w[k];
                    }
                    const auto draw=std::generate_canonical<long double,std::numeric_limits<long double>::digits>(oracle[owners[a]])*total;
                    double cumulative=0;std::size_t chosen=w.size()-1;
                    for(std::size_t k=0;k<w.size();++k) {cumulative+=w[k];if(draw<cumulative) {chosen=k;break;}}
                    expected[a]=chosen%2;if(width==2) expected[b]=chosen/2;
                }
                std::map<int,int> counts;for(auto site:line)counts[site]=2;
                const auto native_plan=tnbp::mcmc_build_edge_line_plan(edges,line,groups,owners,step%2,width,
                    step%3==0?0:step%3==1?2:8,step%3==0?0:step%3==1?8:2,counts);
                tnbp::mcmc_update_line_heat_bath(ctx,edges,W,ids,owners,upper,upper_ids,target,ids,lower,lower_ids,
                    line,groups,state,comm,rng,step%2,width,nullptr,nullptr,&native_plan);
                for(std::size_t k=0;k<ids.size();++k) require(state[k]==expected[ids[k]],"square production sweep differs from roulette oracle");
            }
        }
    }
    if(rank==0) std::cout<<"square D2 chi3 PASS: p2 weighted MAC=612/peak36, p3 common-subtree MAC=756/peak72; plan reuse, all states, widths 1/2, both directions, upper/lower first, distributed owners and empty rank0\n";
}

// Brick-wall representation of a finite honeycomb patch: three four-site
// horizontal chains and alternating upper/lower rungs. The middle targets
// have at most three virtual legs; group members need not touch the target.
void honeycomb_check(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    constexpr int n=4,D=2,chi=3;
    double max_error=0;std::size_t max_work=0,max_elements=0;
    int boundary_first=0,target_first=0;
    for(int parity:{0,1}) {
        Edges edges;std::vector<int> dimensions;
        for(int row=0;row<3;++row) for(int j=0;j+1<n;++j) {
            edges.emplace_back(row*n+j,row*n+j+1);dimensions.push_back(row==0 ? D:chi);
        }
        for(int j=0;j<n;++j) {edges.emplace_back(j,(j%2==parity ? n:2*n)+j);dimensions.push_back(D);}
        // Change both vector order and endpoint orientation without changing
        // topology; fixture tensor axes always follow the resulting edge list.
        if(parity) {
            std::reverse(edges.begin(),edges.end());std::reverse(dimensions.begin(),dimensions.end());
            for(auto& e:edges) std::swap(e.first,e.second);
        }
        const std::vector<int> line{0,1,2,3};
        const std::vector<std::vector<int>> lines{{4,5,6,7},line,{8,9,10,11}};
        const auto groups=tnbp::make_mcmc_grouping(edges,lines).at(1);
        tnbp::validate_mcmc_grouping(edges,lines,tnbp::make_mcmc_grouping(edges,lines));
        bool indirect=false;
        for(const auto& group:groups) for(std::size_t k=1;k<group.size();++k)
            indirect=indirect || tnbp::find_edge(edges,group[0],group[k])==edges.end();
        require(indirect,"honeycomb fixture does not exercise indirect boundary members");
        for(int site:line) require(tnbp::GetSurroundingBondIndex(site,edges).size()<=3,"honeycomb target degree");
        std::size_t combinations=1;for(int d:dimensions) combinations*=d;
        std::vector<Scalar> exact_complex(1<<n,Scalar(0));
        if(rank==0) for(std::size_t config=0;config<combinations;++config) {
            auto index=config;std::vector<std::vector<int>> bonds(3*n);
            for(std::size_t e=0;e<edges.size();++e) {
                const int value=index%dimensions[e];index/=dimensions[e];
                bonds[edges[e].first].push_back(value);bonds[edges[e].second].push_back(value);
            }
            Scalar outer=1.;for(int j=n;j<3*n;++j) outer*=entry(j,0,bonds[j]);
            for(int state=0;state<(1<<n);++state) {
                auto product=outer;for(int j=0;j<n;++j) product*=entry(j,(state>>j)&1,bonds[j]);
                exact_complex[state]+=product;
            }
        }
        std::vector<double> exact_amplitude(1<<n);
        if(rank==0) for(int state=0;state<(1<<n);++state) exact_amplitude[state]=std::abs(exact_complex[state]);
        MPI_Bcast(exact_amplitude.data(),1<<n,MPI_DOUBLE,0,comm);
        for(bool empty_root:{false,true}) {
            std::map<int,int> owners;
            for(int j=0;j<3*n;++j) owners[j]=empty_root ? size-1:(j+1+(size>2 ? j/n:0))%size;
            std::vector<int> ids,upper_ids,lower_ids;std::vector<Tensor> W,upper,lower,target;
            for(int j=0;j<3*n;++j) if(owners[j]==rank) {
                auto bonds=tnbp::GetSurroundingBondIndex(j,edges);
                tcapi::shape_t<Tensor> shape;std::size_t elements=1;
                for(auto e:bonds) {shape.push_back(dimensions[e]);elements*=dimensions[e];}
                if(j<n) shape.push_back(2);
                auto tensor=tcapi::zeros<Tensor>(ctx,shape);
                for(std::size_t v=0;v<elements;++v) for(int state=0;state<(j<n ? 2:1);++state) {
                    auto bits=v;std::vector<int> b;tcapi::elem_coors_t<Tensor> coors;
                    for(auto e:bonds) {b.push_back(bits%dimensions[e]);coors.push_back(b.back());bits/=dimensions[e];}
                    if(j<n) coors.push_back(state);
                    tcapi::set_elem(ctx,tensor,coors,entry(j,state,b));
                }
                if(j<n) {ids.push_back(j);W.push_back(std::move(tensor));shape.pop_back();target.push_back(tcapi::zeros<Tensor>(ctx,shape));}
                else if(j<2*n) {upper_ids.push_back(j);upper.push_back(std::move(tensor));}
                else {lower_ids.push_back(j);lower.push_back(std::move(tensor));}
            }
            std::vector<int> state(ids.size(),0),dims;
            std::vector<std::vector<Tensor>> slices;std::vector<std::vector<int>> axes;
            std::vector<std::vector<std::vector<Tensor>>> boundary;
            std::vector<std::vector<std::vector<int>>> boundary_axes;
            std::vector<std::vector<int>> boundary_owners;
            std::vector<tcapi::bond_dim_t<Tensor>> bond_dims;
            tnbp::mcmc_prepare_heat_bath_cache(ctx,edges,W,ids,owners,upper,upper_ids,lower,lower_ids,
                line,groups,state,comm,slices,axes,dims,bond_dims,nullptr,&boundary,&boundary_axes,&boundary_owners);
            std::vector<int> group_owners;for(int j:line) group_owners.push_back(owners[j]);
            for(int initial=0;initial<(1<<n);++initial) for(int direction:{0,1}) for(int width:{1,2}) {
                std::vector<int> selected(n);for(int j:line) selected[j]=(initial>>j)&1;
                std::vector<std::size_t> order{0,1,2,3};if(direction) std::reverse(order.begin(),order.end());
                std::vector<std::vector<Tensor>> suffix;std::vector<std::vector<int>> suffix_axes;
                std::vector<std::size_t> audit(10,0);
                tnbp::mcmc_build_heat_bath_environments(ctx,slices,axes,group_owners,selected,order,bond_dims,
                    comm,suffix,suffix_axes,width,&audit,&boundary,&boundary_axes,&boundary_owners);
                // Retain full cache fingerprints across window evaluations.
                auto fingerprint=[&]() {
                    std::vector<Scalar> values;
                    for(const auto& cache:suffix) for(const auto& t:cache) {
                        auto shape=tcapi::shape(ctx,t);std::size_t count=1;for(auto d:shape) count*=d;
                        for(std::size_t v=0;v<count;++v) {
                            auto index=v;tcapi::elem_coors_t<Tensor> coors;
                            for(auto d:shape) {coors.push_back(index%d);index/=d;}
                            values.push_back(tcapi::get_elem(ctx,t,coors));
                        }
                    }
                    return values;
                };
                const auto before=fingerprint();const auto environment_calls=audit[0];
                std::vector<Tensor> prefix;std::vector<int> prefix_axes;
                for(std::size_t pos=0;pos+width<=n;++pos) {
                    const auto a=order[pos],b=width==2 ? order[pos+1]:a,end=pos+width;
                    // Trace an actual production candidate-helper call using
                    // the proper advancing side (before changing that side).
                    if(initial==0 && width==1 && pos>0) {
                        std::vector<std::vector<Tensor>> traced;std::vector<int> traced_axes,steps;
                        std::vector<std::vector<std::size_t>> trace;
                        tnbp::mcmc_prepare_group_candidates(ctx,slices[a],axes[a],owners[a],2,
                            boundary[a],boundary_axes[a],boundary_owners[a],prefix,prefix_axes,true,owners[order[pos-1]],
                            bond_dims,comm,traced,traced_axes,1,nullptr,&steps,0,&trace);
                        require(!steps.empty(),"honeycomb candidate trace missing");
                        if(steps.front()==-1) ++target_first;else ++boundary_first;
                        for(const auto& op:trace) {max_work=std::max(max_work,op[3]);max_elements=std::max(max_elements,op[4]);}
                        for(const auto& t:traced) require(t.size()==std::size_t(rank==owners[a]),"honeycomb trace final owner");
                    }
                    std::vector<double> amplitudes,weights;std::vector<std::vector<Tensor>> candidates;std::vector<int> candidate_axes;
                    tnbp::mcmc_evaluate_window_weights(ctx,slices[a],axes[a],owners[a],2,
                        slices[b],axes[b],owners[b],width==2 ? 2:1,width==2,
                        prefix,prefix_axes,pos>0,pos>0 ? owners[order[pos-1]]:owners[a],
                        suffix[end],suffix_axes[end],end<n,end<n ? owners[order[end]]:owners[a],
                        bond_dims,comm,amplitudes,weights,candidates,candidate_axes,&audit,
                        boundary[a],boundary_axes[a],boundary_owners[a],boundary[b],boundary_axes[b],boundary_owners[b]);
                    double reference_total=0,weight_total=std::accumulate(weights.begin(),weights.end(),0.);
                    std::vector<double> reference(amplitudes.size());
                    for(std::size_t k=0;k<amplitudes.size();++k) {
                        auto c=selected;c[a]=k%2;if(width==2) c[b]=k/2;
                        int encoded=0;for(int j:line) encoded+=c[j]<<j;
                        const auto expected=exact_amplitude[encoded];
                        max_error=std::max(max_error,std::abs(amplitudes[k]-expected)/(1+expected));
                        require(std::abs(amplitudes[k]-expected)<1.e-11*(1+expected),"honeycomb candidate amplitude mismatch");
                        reference[k]=expected*expected;reference_total+=reference[k];
                    }
                    for(std::size_t k=0;k<reference.size();++k)
                        require(std::abs(reference[k]/reference_total-weights[k]/weight_total)<1.e-11,"honeycomb conditional probability mismatch");
                    require(audit[0]==environment_calls && audit[5]==1 && audit[6]==0,"honeycomb rebuilt exterior or used identities");
                    for(const auto& t:candidates) require(t.size()==std::size_t(rank==owners[a]),"honeycomb candidate owner");
                    const auto* storage=candidates[selected[a]].data();
                    tnbp::mcmc_absorb_group(std::size_t(selected[a]),candidates,candidate_axes,prefix,prefix_axes);
                    require(prefix.data()==storage,"honeycomb selected intermediate not moved");
                }
                require(fingerprint()==before,"honeycomb candidate mutated exterior cache");
            }
            for(int width:{1,2}) {
                std::fill(state.begin(),state.end(),0);std::vector<int> expected(n,0);
                std::mt19937 rng(631+rank);std::vector<std::mt19937> oracle;
                for(int r=0;r<size;++r) oracle.emplace_back(631+r);
                for(int step=0;step<16;++step) {
                    auto order=line;if(step%2) std::reverse(order.begin(),order.end());
                    for(int pos=0;pos+width<=n;++pos) {
                        const int a=order[pos],b=width==2 ? order[pos+1]:a;
                        std::vector<double> weights(width==2 ? 4:2);double total=0;
                        for(std::size_t k=0;k<weights.size();++k) {
                            auto c=expected;c[a]=k%2;if(width==2) c[b]=k/2;
                            int encoded=0;for(int j:line) encoded+=c[j]<<j;
                            const auto amp=exact_amplitude[encoded];weights[k]=amp*amp;total+=weights[k];
                        }
                        const auto draw=std::generate_canonical<long double,std::numeric_limits<long double>::digits>(oracle[owners[a]])*total;
                        double cumulative=0;std::size_t chosen=weights.size()-1;
                        for(std::size_t k=0;k<weights.size();++k) {cumulative+=weights[k];if(draw<cumulative) {chosen=k;break;}}
                        expected[a]=chosen%2;if(width==2) expected[b]=chosen/2;
                    }
                    std::map<int,int> counts;for(auto site:line)counts[site]=2;
                    const auto native_plan=tnbp::mcmc_build_edge_line_plan(edges,line,groups,owners,step%2,width,
                        step%3==0?0:step%3==1?2:8,step%3==0?0:step%3==1?8:2,counts);
                    tnbp::mcmc_update_line_heat_bath(ctx,edges,W,ids,owners,upper,upper_ids,target,ids,lower,lower_ids,
                        line,groups,state,comm,rng,step%2,width,nullptr,nullptr,&native_plan);
                    for(std::size_t k=0;k<ids.size();++k) require(state[k]==expected[ids[k]],"honeycomb sweep differs from independent roulette");
                }
            }
        }
    }
    require(boundary_first+target_first>0,"honeycomb plan trace missing");
    if(rank==0) std::cout<<"honeycomb D2 chi3 PASS: two hexagons, 16 configurations, all candidates/windows/directions, both rung parities, planned indirect members, empty root, immutable caches/move/roulette; max relative amplitude error "
        <<max_error<<", traced max MAC "<<max_work<<", max intermediate elements "<<max_elements
        <<", boundary-first "<<boundary_first<<", target-first "<<target_first<<"\n";
}

// Unit-cell graph A(m,n)--B(m,n), A(m,n)--B(m-1,n),
// A(m,n)--B(m,n-1). Test actual BFS layers and their paired zigzag chains.
void zigzag_check(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    Edges edges;
    for(int row=0;row<2;++row) for(int cell=0;cell<3;++cell) {
        const int a=6*row+2*cell;edges.emplace_back(a,a+1);
        if(cell) edges.emplace_back(a,a-1);
        if(row) edges.emplace_back(a,a-5);
    }
    const std::vector<int> initial_boundary{0,2,4};
    std::vector<std::vector<int>> automatic;
    tnbp::bmps_setup_lines(edges,initial_boundary,automatic);
    require(automatic==std::vector<std::vector<int>>{{0,2,4},{1,3,5},{6,8,10},{7,9,11}},"unexpected automatic zigzag layers");
    const std::vector<std::vector<int>> paired{{0,1,2,3,4,5},{6,7,8,9,10,11}};
    std::vector<int> dims(12,2),free;
    for(int site:{0,2,4,7,9,11}) dims[site]=1;
    for(int site=0;site<12;++site) if(dims[site]==2) free.push_back(site);
    std::vector<double> reference(1<<12,0.);
    if(rank==0) for(int c=0;c<(1<<free.size());++c) {
        std::vector<int> state(12,0);int encoded=0;
        for(std::size_t k=0;k<free.size();++k) {state[free[k]]=(c>>k)&1;encoded|=state[free[k]]<<free[k];}
        reference[encoded]=std::abs(exact(edges,state));
    }
    MPI_Bcast(reference.data(),reference.size(),MPI_DOUBLE,0,comm);
    double max_error=0.;std::size_t checked=0;
    for(bool empty_root:{false,true}) {
        std::map<int,int> owners;std::vector<int> ids;
        for(int site=0;site<12;++site) {
            owners[site]=empty_root ? size-1:(site+1)%size;
            if(owners[site]==rank) ids.push_back(site);
        }
        const auto W=tensors(ctx,edges,dims,ids);
        for(bool pair_layers:{false,true}) {
            const auto& lines=pair_layers ? paired:automatic;
            const auto plan=tnbp::make_mcmc_grouping(edges,lines);
            for(int c=0;c<(1<<free.size());++c) {
                int encoded=0;for(std::size_t k=0;k<free.size();++k) encoded|=((c>>k)&1)<<free[k];
                std::vector<int> initial;for(int site:ids) initial.push_back((encoded>>site)&1);
                tnbp::MCMCWorkspace<Tensor,int> work;work.lines=lines;
                // Also exercise explicit plan override without changing topology.
                if(c%2) work.manual_grouping=plan;
                tnbp::mcmc_initialize(ctx,edges,W,ids,owners,std::vector<std::vector<int>>{initial},
                    30,1.e-12,tcapi::bond_dim_t<Tensor>(16),1.e-14,1.e-14,work,comm);
                require(work.site_line_sets==plan,"zigzag initialization changed grouping");
                const std::vector<Tensor> empty;const std::vector<int> ie;
                for(std::size_t l=0;l<lines.size();++l) {
                    const auto& line=lines[l];const auto& groups=work.site_line_sets[l];
                    const auto& upper=l+1<lines.size() ? work.T_upper[0][l+1]:empty;
                    const auto& upper_ids=l+1<lines.size() ? work.site_idx_upper[l+1]:ie;
                    const auto& lower=l ? work.T_lower[0][l-1]:empty;
                    const auto& lower_ids=l ? work.site_idx_lower[l-1]:ie;
                    std::vector<std::vector<Tensor>> slices;std::vector<std::vector<int>> axes;
                    std::vector<std::vector<std::vector<Tensor>>> boundary;
                    std::vector<std::vector<std::vector<int>>> boundary_axes;
                    std::vector<std::vector<int>> boundary_owners;
                    std::vector<int> physical,group_owners,selected;
                    for(int site:line) {group_owners.push_back(owners[site]);selected.push_back((encoded>>site)&1);}
                    std::vector<tcapi::bond_dim_t<Tensor>> bond_dims;
                    tnbp::mcmc_prepare_heat_bath_cache(ctx,work.line_edges[l],W,ids,owners,upper,upper_ids,lower,lower_ids,
                        line,groups,initial,comm,slices,axes,physical,bond_dims,nullptr,&boundary,&boundary_axes,&boundary_owners);
                    // BMPS truncation normalizes its tensors. One common
                    // exterior scale per line may differ from the raw network;
                    // it must be identical for every candidate/window/direction.
                    double boundary_scale=-1.;
                    for(int direction:{0,1}) for(int width:{1,2}) {
                        std::vector<std::size_t> order(line.size());std::iota(order.begin(),order.end(),0);
                        if(direction) std::reverse(order.begin(),order.end());
                        std::vector<std::vector<Tensor>> suffix;std::vector<std::vector<int>> suffix_axes;
                        tnbp::mcmc_build_heat_bath_environments(ctx,slices,axes,group_owners,selected,order,bond_dims,
                            comm,suffix,suffix_axes,width,nullptr,&boundary,&boundary_axes,&boundary_owners);
                        std::vector<Tensor> prefix;std::vector<int> prefix_axes;
                        for(std::size_t pos=0;pos+width<=line.size();++pos) {
                            const auto a=order[pos],b=width==2 ? order[pos+1]:a,end=pos+width;
                            std::vector<double> amplitudes,weights;std::vector<std::vector<Tensor>> candidates;std::vector<int> candidate_axes;
                            tnbp::mcmc_evaluate_window_weights(ctx,slices[a],axes[a],owners[line[a]],physical[a],
                                slices[b],axes[b],owners[line[b]],width==2 ? physical[b]:1,width==2,
                                prefix,prefix_axes,pos>0,pos>0 ? owners[line[order[pos-1]]]:owners[line[a]],
                                suffix[end],suffix_axes[end],end<line.size(),end<line.size() ? owners[line[order[end]]]:owners[line[a]],
                                bond_dims,comm,amplitudes,weights,candidates,candidate_axes,nullptr,
                                boundary[a],boundary_axes[a],boundary_owners[a],boundary[b],boundary_axes[b],boundary_owners[b]);
                            for(std::size_t k=0;k<amplitudes.size();++k) {
                                int candidate=(encoded & ~(1<<line[a])) | (int(k%physical[a])<<line[a]);
                                if(width==2) candidate=(candidate & ~(1<<line[b])) | (int(k/physical[a])<<line[b]);
                                const auto expected=reference[candidate];
                                if(boundary_scale<0 && expected>0) boundary_scale=amplitudes[k]/expected;
                                require(std::isfinite(boundary_scale) && boundary_scale>0,"invalid common BMPS scale");
                                const auto scaled=amplitudes[k]/boundary_scale;
                                max_error=std::max(max_error,std::abs(scaled-expected)/(1+expected));
                                require(std::abs(scaled-expected)<1.e-10*(1+expected),"zigzag BMPS candidate differs from original graph up to common exterior scale");
                                ++checked;
                            }
                            const auto* storage=candidates[selected[a]].data();
                            tnbp::mcmc_absorb_group(std::size_t(selected[a]),candidates,candidate_axes,prefix,prefix_axes);
                            require(prefix.data()==storage,"zigzag selected storage not moved");
                        }
                        // Independent whole-line roulette oracle, with the
                        // exterior fixed to this exact initial configuration.
                        int expected=encoded;std::vector<std::mt19937> oracle;
                        for(int r=0;r<size;++r) oracle.emplace_back(811+r);
                        for(std::size_t pos=0;pos+width<=line.size();++pos) {
                            const auto a=order[pos],b=width==2 ? order[pos+1]:a;
                            std::vector<double> weights(physical[a]*(width==2 ? physical[b]:1));double total=0;
                            for(std::size_t k=0;k<weights.size();++k) {
                                auto candidate=(expected & ~(1<<line[a])) | (int(k%physical[a])<<line[a]);
                                if(width==2) candidate=(candidate & ~(1<<line[b])) | (int(k/physical[a])<<line[b]);
                                weights[k]=reference[candidate]*reference[candidate];total+=weights[k];
                            }
                            const auto draw=std::generate_canonical<long double,std::numeric_limits<long double>::digits>(oracle[owners[line[a]]])*total;
                            double cumulative=0;std::size_t chosen=weights.size()-1;
                            for(std::size_t k=0;k<weights.size();++k) {cumulative+=weights[k];if(draw<cumulative) {chosen=k;break;}}
                            expected=(expected & ~(1<<line[a])) | (int(chosen%physical[a])<<line[a]);
                            if(width==2) expected=(expected & ~(1<<line[b])) | (int(chosen/physical[a])<<line[b]);
                        }
                        auto states=initial;std::vector<Tensor> target;std::vector<int> target_ids;
                        for(std::size_t k=0;k<ids.size();++k) if(std::find(line.begin(),line.end(),ids[k])!=line.end()) {
                            target_ids.push_back(ids[k]);target.push_back(tcapi::copy(ctx,work.V[0][k]));
                        }
                        std::mt19937 rng(811+rank);
                        std::map<int,int> counts;for(auto site:line)counts[site]=dims[site];
                        const auto native_plan=tnbp::mcmc_build_edge_line_plan(work.line_edges[l],line,groups,owners,direction,width,2,16,counts,
                            tnbp::mcmc_boundary_edge_model(work.line_edges[l],lines,l));
                        tnbp::mcmc_update_line_heat_bath(ctx,work.line_edges[l],W,ids,owners,upper,upper_ids,target,target_ids,lower,lower_ids,
                            line,groups,states,comm,rng,direction,width,nullptr,nullptr,&native_plan);
                        for(std::size_t k=0;k<ids.size();++k) require(states[k]==((expected>>ids[k])&1),"zigzag whole-line update differs from roulette oracle");
                    }
                }
            }
        }
    }
    if(rank==0) std::cout<<"zigzag unit-cell honeycomb PASS: auto A/B layers and two paired zigzag lines; 64 physical configurations, "
        <<checked<<" candidate amplitudes, all lines/windows/widths/directions, actual BMPS, manual/auto grouping, owner roulette, empty root; max relative amplitude error "<<max_error<<"\n";
}

void native_scalar_metadata_check(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    std::map<int,int> owners{{0,size-1},{1,0}};std::vector<int> members;
    for(int r=0;r<size;++r)members.push_back(r);
    tnbp::mcmc_edge_region_plan<int> plan{Edges{{0,1}},{{0,1,2}},{{0,1,0}},members,0,{},0,{}, {}};
    tnbp::mcmc_compile_edge_region(plan);
    std::vector<Tensor> a,b;
    for(int site:{0,1})if(rank==owners[site])for(int state=0;state<(site==0?3:2);++state) {
        auto value=tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{2,3,4});
        for(int i=0;i<2;++i)for(int j=0;j<3;++j)for(int k=0;k<4;++k)
            tcapi::set_elem(ctx,value,tcapi::elem_coors_t<Tensor>{i,j,k},entry(site,state,{i,j,k}));
        (site==0?a:b).push_back(std::move(value));
    }
    const auto reductions=allreduce_calls,broadcasts=bcast_calls;
    const auto metadata=tnbp::mcmc_prepare_edge_scalar_join(ctx,plan,owners,members,rank,
        rank==owners[0]?&a.at(0):nullptr,rank==owners[1]?&b.at(0):nullptr);
    const auto saved=metadata;std::vector<std::size_t> audit(15,0);
    for(int t=0;t<2;++t)for(int s=0;s<3;++s) {
        auto out=tnbp::mcmc_execute_edge_scalar_join(ctx,plan,metadata,
            rank==owners[0]?&a.at(s):nullptr,rank==owners[1]?&b.at(t):nullptr,s==2,comm,&audit);
        require(out.size()==std::size_t(rank==owners[0]),"scalar result owner");
        if(rank==owners[0]) {
            Scalar expected=0;
            for(int i=0;i<2;++i)for(int j=0;j<3;++j)for(int k=0;k<4;++k)
                expected+=entry(0,s,{i,j,k})*entry(1,t,{i,j,k});
            require(std::abs(tcapi::get_elem(ctx,out[0],tcapi::elem_coors_t<Tensor>{})-expected)<1.e-11*(1+std::abs(expected)),"prepared scalar join changed composite-factor order");
            require(tcapi::shape(ctx,a[s])==tcapi::shape_t<Tensor>({2,3,4}),"roulette-pending A was consumed or packed in place");
        }
    }
    require(metadata==saved,"scalar descriptor changed across candidates");
    require(allreduce_calls==reductions && bcast_calls==broadcasts,"scalar descriptor/execution added collective");
    require(audit[11]==std::size_t(rank==owners[1]?2:0),"RHS last-use moves changed");
    require(audit[7]==std::size_t((rank==owners[0]?6:0)+(rank==owners[1]?4:0)),"scalar input copy count changed");
    if(rank==0)std::cout<<"prepared scalar metadata PASS: six unequal-factor complex joins, descriptor reuse, A protected, RHS last-use move, zero collectives\n";
}

void boundary_model_check() {
    const std::vector<std::vector<int>> lines{{0,1},{2,3},{4,5},{6,7}};
    Edges edges{{0,1},{2,3},{4,5},{6,7},{0,2},{1,3},{2,4},{3,5},{4,6},{5,7}};
    for(std::size_t k=0;k<lines.size();++k) {
        Edges reduced;tnbp::mcmc_init_edges(edges,lines,static_cast<int>(k),reduced);
        const auto model=tnbp::mcmc_boundary_edge_model(reduced,lines,k);
        const auto truncated=k%2==0?std::pair<int,int>{2,3}:std::pair<int,int>{4,5};
        bool found=false;
        for(std::size_t e=0;e<reduced.size();++e) {
            const bool boundary=tnbp::make_edge(reduced[e].first,reduced[e].second)==truncated;
            require(model[e]==(boundary?"x":"D"),"boundary checkpoint confused raw D with truncated chi");
            found=found || boundary;
        }
        require(found,"boundary model fixture missing truncated bond");
    }
}

// A common subtree fuses both the target's axes and an outside neighbor's
// axes. Unequal factors and interleaved edges expose reversed fusion order.
void native_region_check(Context& ctx,MPI_Comm comm) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    const Edges edges{{1,3},{0,2},{0,1},{2,3},{1,2},{0,3}};
    const std::vector<int> dims{3,2,2,4,3,2};
    std::map<int,int> owners;std::vector<int> members;
    for(int site=0;site<4;++site)owners[site]=size==3?(site==2?2:site==3?1:0):site%size;
    for(int r=0;r<size;++r)members.push_back(r);
    tnbp::mcmc_edge_ports ports;for(std::size_t e=0;e<edges.size();++e)ports.push_back({e});
    tnbp::mcmc_edge_region_plan<int> plan{edges,ports,{{0,1,0},{2,0,2}},members,2,{},0,{}, {}};
    tnbp::mcmc_compile_edge_region(plan);
    std::map<int,std::vector<Tensor>> raw;
    for(int site=0;site<4;++site)if(owners[site]==rank) {
        const auto bonds=tnbp::GetSurroundingBondIndex(site,edges);
        tcapi::shape_t<Tensor> shape;std::size_t n=1;
        for(auto e:bonds){shape.push_back(dims[e]);n*=dims[e];}
        for(int s=0;s<(site==2?3:1);++s) {
            auto tensor=tcapi::zeros<Tensor>(ctx,shape);
            for(std::size_t i=0;i<n;++i) {
                auto bits=i;tcapi::elem_coors_t<Tensor> coors;std::vector<int> values;
                for(auto d:shape){coors.push_back(bits%d);values.push_back(bits%d);bits/=d;}
                tcapi::set_elem(ctx,tensor,coors,entry(site,s,values));
            }
            raw[site].push_back(std::move(tensor));
        }
    }
    const auto reductions=allreduce_calls,broadcasts=bcast_calls;
    std::vector<std::size_t> audit(15,0);
    auto provider=[&](int site,int state)->const Tensor&{return raw.at(site).at(site==2?state:0);};
    auto values=tnbp::mcmc_execute_edge_region<Tensor>(ctx,plan,owners,members,comm,3,provider,1,&audit);
    require(allreduce_calls==reductions && bcast_calls==broadcasts,"native graph executor used a collective");
    require(audit[8]==1 && audit[1]==4,"common subtree was not reused across target states");
    unsigned long long copies=audit[7],total_copies=0;
    MPI_Allreduce(&copies,&total_copies,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,comm);
    require(total_copies==10,"native branch inputs copied again between steps, or common values not moved at last use");
    require(values.size()==std::size_t(rank==owners[2]?3:0),"native root changed target owner");
    if(rank==owners[2])for(int s=0;s<3;++s)for(int x=0;x<3;++x)for(int y=0;y<4;++y)for(int z=0;z<2;++z) {
        Scalar expected=0;
        for(int a=0;a<2;++a)for(int b=0;b<2;++b)for(int c=0;c<3;++c)
            expected+=entry(0,0,{a,b,z})*entry(1,0,{x,b,c})*entry(2,s,{a,y,c});
        const auto actual=tcapi::get_elem(ctx,values[s],tcapi::elem_coors_t<Tensor>{x,y,z});
        require(std::abs(actual-expected)<1.e-11*(1+std::abs(expected)),"native common-neighbor/cut factor order");
    }
    // Raw target and exterior neighbor are immutable despite native fusions.
    for(const auto& item:raw)for(std::size_t state=0;state<item.second.size();++state) {
        const auto sh=tcapi::shape(ctx,item.second[state]);
        const auto bonds=tnbp::GetSurroundingBondIndex(item.first,edges);
        require(sh.size()==bonds.size(),"native source cache was fused in place");
        for(std::size_t i=0;i<sh.size();++i)require(sh[i]==dims[bonds[i]],"native cache factor dimension changed");
        require(tcapi::get_elem(ctx,item.second[state],tcapi::elem_coors_t<Tensor>(sh.size(),0))==entry(item.first,state,std::vector<int>(sh.size(),0)),"native cache value changed");
    }
    if(rank==0)std::cout<<"native original-edge region PASS: outside neighbor, unequal/interleaved fusion, all three states, common subtree once, owner/cut/cache invariant, zero collectives\n";
}

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int result=0;
    try {
        initialize_tensor_device(MPI_COMM_WORLD);Context ctx;tcapi::create_context(ctx);
        if(argc>2 && std::string(argv[1])=="fatal") {
            fatal_failure_check(ctx,MPI_COMM_WORLD,std::stoi(argv[2]));
            tcapi::destroy_context(ctx);MPI_Finalize();return 0;
        }
        commit_move_only_check(MPI_COMM_WORLD);
        if(argc>1 && std::string(argv[1])=="zigzag") {
            zigzag_check(ctx,MPI_COMM_WORLD);tcapi::destroy_context(ctx);MPI_Finalize();return 0;
        }
        if(argc>1 && std::string(argv[1])=="honeycomb") {
            honeycomb_check(ctx,MPI_COMM_WORLD);tcapi::destroy_context(ctx);MPI_Finalize();return 0;
        }
        Edges k4{{0,1},{2,0},{3,0},{1,2},{3,1},{2,3}};
        if(rank==0) roulette_check(ctx);
        conditional_check(ctx,k4,{2,2,2,2},MPI_COMM_WORLD);
        conditional_check(ctx,k4,{2,3,2,1},MPI_COMM_WORLD);
        distributed_check(ctx,k4,{2,3,2,1},MPI_COMM_WORLD);
        support_check(ctx,MPI_COMM_WORLD);
        native_scalar_metadata_check(ctx,MPI_COMM_WORLD);
        boundary_model_check();
        native_region_check(ctx,MPI_COMM_WORLD);
        cache_reuse_check(ctx,MPI_COMM_WORLD);
        line_communicator_check(ctx,MPI_COMM_WORLD);
        review_fixes_check(ctx,MPI_COMM_WORLD);
        planner_fusion_check(ctx,MPI_COMM_WORLD);
        scalar_failure_check(ctx,MPI_COMM_WORLD);
        square_order_check(ctx,MPI_COMM_WORLD);
        tcapi::destroy_context(ctx);
        if(rank==0) std::cout<<"full-state heat-bath regression PASS\n";
    }catch(const std::exception& e) {std::cerr<<"rank "<<rank<<": "<<e.what()<<"\n";MPI_Abort(MPI_COMM_WORLD,1);result=1;}
    MPI_Finalize();return result;
}
