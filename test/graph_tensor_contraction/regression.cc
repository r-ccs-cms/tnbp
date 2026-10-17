#include "../common/device.h"
#include "../common/graph_fixture.h"
template<class S>
void run(tcapi::context_handle_t<TensorOf<S>>& ctx,const std::string& name,Edges edges,const std::vector<int>& sites,
         int ownership_shift,bool zero,MPI_Comm comm,bool single) {
    run_shift_fixture<S>(ctx,name,edges,sites,ownership_shift,zero,comm,
        [](const auto&,auto&){}, [](auto&...){}, single);
}
// Explicit edge-order expectations, independent of the shared plan used by
// both public paths. The B edge precedes the retained A edge for neighbor 2.
void graph_plan_checks() {
    const Edges original{{1,2},{7,9},{0,4},{0,2},{1,4},{1,5},{1,0},{0,6}};
    for(bool reverse:{false,true}) {
        auto edges=original;
        Edges expected{{7,9},{4,8},{2,8},{5,8},{6,8}};
        if(reverse) {std::reverse(edges.begin(),edges.end());std::reverse(expected.begin(),expected.end());}
        const auto before=edges;
        const auto plan=tnbp::detail::make_graph_contraction_plan(edges,0,1,8);
        if(edges!=before || plan.edges!=expected || plan.fusions.size()!=2)
            throw std::runtime_error("graph contraction plan changes edge order");
        tnbp::graph_tensor_contraction_for_edges(edges,0,1,8);
        if(edges!=expected) throw std::runtime_error("graph-only edge order mismatch");
    }
    Edges scalar{{1,0}};
    tnbp::graph_tensor_contraction_for_edges(scalar,0,1,0);
    if(!scalar.empty()) throw std::runtime_error("scalar contraction retains an edge");
    for(int kind=0;kind<4;++kind) {
        auto edges=original;int a=0,b=1,c=8;
        if(kind==0) b=7;
        if(kind==1) c=9;
        if(kind==2) edges.push_back({2,1});
        if(kind==3) c=1;
        const auto before=edges;bool threw=false;
        try {tnbp::graph_tensor_contraction_for_edges(edges,a,b,c);}
        catch(const std::invalid_argument&) {threw=true;}
        if(!threw || edges!=before) throw std::runtime_error("invalid graph-only plan mutated edges");
    }
}

// Direct contraction coverage: two common neighbors, unequal bond dimensions,
// both edge orders, partial tensors/owners and fresh/result-in-place labels.
template<class S>
void contraction_checks(MPI_Comm comm) {
    using T=TensorOf<S>;tcapi::context_handle_t<T> ctx;tcapi::create_context(ctx);
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    for(bool reverse:{false,true}) for(int result:{0,8}) {
        Edges edges{{0,1},{0,2},{1,3},{1,2},{0,3},{0,10},{1,11},{2,12},{3,13}};
        if(reverse) std::reverse(edges.begin(),edges.end());
        std::vector<int> ids;std::vector<T> ts;std::map<int,int> owners;
        for(int site=0;site<4;++site) {
            owners[site]=site%size;if(site%size!=rank) continue;
            ids.push_back(site);tcapi::shape_t<T> shape;
            for(auto e:tnbp::GetSurroundingBondIndex(site,edges)) shape.push_back(e%3==0 ? 3 : 2);
            auto t=tcapi::zeros<T>(ctx,shape);
            std::vector<S> data(tcapi::size(ctx,t));
            for(std::size_t j=0;j<data.size();++j) data[j]=value<S>(site,j);
            auto it=data.begin();
            t=tcapi::assign_from_range<T>(ctx,shape,it,[shape](const auto& c){return tnbp::address_from_coor(shape,c);});
            ts.push_back(std::move(t));
        }
        const auto expected=boundary(ctx,{0,1,2,3},edges,ts,ids,owners,comm);
        auto planned=edges;tnbp::graph_tensor_contraction_for_edges(planned,0,1,result);
        tnbp::graph_tensor_contraction(ctx,edges,ts,ids,owners,0,1,result,comm);
        if(edges!=planned || owners.count(1) || owners.at(result)!=0)
            throw std::runtime_error("contraction graph/owner mismatch");
        const auto actual=boundary(ctx,{result,2,3},edges,ts,ids,owners,comm);
        if(actual.size()!=expected.size()) throw std::runtime_error("contraction output shape mismatch");
        for(std::size_t j=0;j<actual.size();++j)
            if(std::abs(actual[j]-expected[j])>1.e-10*(1+std::abs(expected[j])))
                throw std::runtime_error("contraction open-boundary value mismatch");
    }
    if(rank==0) std::cout<<"direct contraction fusion/edge order/partial owners PASS"<<std::endl;
}

void invalid_contraction_checks(MPI_Comm comm) {
    using T=TensorOf<double>;tcapi::context_handle_t<T> ctx;tcapi::create_context(ctx);
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    for(int kind=0;kind<15;++kind) {
        Edges edges{{0,1},{0,2},{1,2},{0,10}};
        std::map<int,int> owners{{0,0},{1,1%size},{2,2%size}};
        std::vector<int> ids;std::vector<T> ts;
        for(int site=0;site<3;++site) if(owners.at(site)==rank) {
            ids.push_back(site);
            const auto degree=tnbp::GetSurroundingBondIndex(site,edges).size();
            auto t=tcapi::zeros<T>(ctx,tcapi::shape_t<T>(degree,2));
            tcapi::set_elem(ctx,t,tcapi::elem_coors_t<T>(degree,0),3.);
            ts.push_back(std::move(t));
        }
        int a=0,b=1,c=0;
        if(kind==0) b=a;
        if(kind==1) b=10; // nonadjacent to 1, so choose a=1 below
        if(kind==1) {a=1;c=1;}
        if(kind==2) c=1; // result cannot overwrite b
        if(kind==3) c=10; // external graph label, absent from owners
        if(kind==4) edges.push_back({1,0});
        if(kind==5) edges.push_back({2,2});
        if(kind==6) owners.erase(2); // missing common-neighbor owner
        if(kind==7) owners[2]=size; // invalid rank
        if(kind==8 && rank==2%size) {
            auto i=std::find(ids.begin(),ids.end(),2)-ids.begin();
            ids.erase(ids.begin()+i);ts.erase(ts.begin()+i);
        }
        if(kind==9 && rank==0) ts[0]=tcapi::zeros<T>(ctx,{2}); // wrong axis count on one rank
        if(kind==10 && rank==1%size) {
            auto i=std::find(ids.begin(),ids.end(),1)-ids.begin();
            ts[i]=tcapi::zeros<T>(ctx,{3,2}); // contraction bond disagrees across ranks
        }
        if(kind==11 && rank==0) ids.push_back(0); // local storage mismatch
        if(kind==12) owners[0]=-1; // public API does not silently skip operands
        if(kind==13) owners[2]=-1; // or missing affected boundary tensors
        if(kind==14) {owners[8]=0;c=8;} // isolated owned label collision
        const auto original_edges=edges;const auto original_ids=ids;const auto original_owners=owners;
        std::vector<tcapi::shape_t<T>> shapes;std::vector<std::vector<double>> values;
        for(const auto& t:ts) {
            auto shape=tcapi::shape(ctx,t);shapes.push_back(shape);
            values.emplace_back(tcapi::size(ctx,t));
            tcapi::to_range(ctx,t,values.back().begin(),[shape](const auto& co){return tnbp::address_from_coor(shape,co);});
        }
        bool threw=false;
        try {tnbp::graph_tensor_contraction(ctx,edges,ts,ids,owners,a,b,c,comm);}
        catch(const std::invalid_argument&) {threw=true;}
        if(!threw || edges!=original_edges || ids!=original_ids || owners!=original_owners || ts.size()!=values.size())
            throw std::runtime_error("contraction preflight failed to reject before mutation");
        for(std::size_t i=0;i<ts.size();++i) {
            const auto shape=tcapi::shape(ctx,ts[i]);std::vector<double> data(tcapi::size(ctx,ts[i]));
            tcapi::to_range(ctx,ts[i],data.begin(),[shape](const auto& co){return tnbp::address_from_coor(shape,co);});
            if(shape!=shapes[i] || data!=values[i]) throw std::runtime_error("invalid contraction mutated tensors");
        }
    }
    if(rank==0) std::cout<<"collective contraction preflight: 15 invalid cases PASS"<<std::endl;
}

void invalid_checks(MPI_Comm comm) {
    using T=TensorOf<double>;tcapi::context_handle_t<T> ctx;tcapi::create_context(ctx);
    int rank;MPI_Comm_rank(comm,&rank);
    for(int kind=0;kind<5;++kind) {
        Edges edges{{0,2},{1,2}};std::vector<T> ts;std::vector<int> ids;
        std::map<int,int> owners;
        bool threw=false;
        try {
            if(kind==0) tnbp::graph_shift_bond_to_neighbor(ctx,edges,ts,ids,owners,0,0,2,comm);
            if(kind==1) tnbp::graph_shift_bond_to_neighbor(ctx,edges,ts,ids,owners,0,2,2,comm);
            if(kind==2) tnbp::graph_shift_bond_to_neighbor(ctx,edges,ts,ids,owners,0,1,99,comm);
            if(kind==3) tnbp::graph_shift_bonds_to_neighbor(ctx,edges,ts,ids,owners,0,1,std::vector<int>{},comm);
            if(kind==4) tnbp::graph_shift_bonds_to_neighbor(ctx,edges,ts,ids,owners,0,1,std::vector<int>{2,2},comm);
        } catch(const std::invalid_argument&) {threw=true;}
        if(!threw || edges!=Edges{{0,2},{1,2}}) throw std::runtime_error("invalid selected bond accepted");
    }
}

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    initialize_tensor_device(MPI_COMM_WORLD);
    try {
        graph_plan_checks();
        invalid_checks(MPI_COMM_WORLD);
        invalid_contraction_checks(MPI_COMM_WORLD);
        auto cases=[&](auto scalar) {
            using S=decltype(scalar);
            tcapi::context_handle_t<TensorOf<S>> ctx;tcapi::create_context(ctx);
            contraction_checks<S>(MPI_COMM_WORLD);
            for(int shift:{0,1}) {
                run<S>(ctx,"single partial, no neighbor",{{0,3},{1,3},{0,10},{1,11},{2,12},{3,13}},
                       {0,1,2,3},shift,false,MPI_COMM_WORLD,true);
                run<S>(ctx,"single partial, shared endpoint",{{0,3},{1,3},{2,3},{0,2},{0,10},{1,11},{2,12},{3,13}},
                       {0,1,2,3},shift,false,MPI_COMM_WORLD,true);
                run<S>(ctx,"single scalar remainder",{{0,2}}, {0,1,2},shift,false,MPI_COMM_WORLD,true);
            }
        };
        cases(double{});cases(std::complex<double>{});
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,1);
    }
    MPI_Finalize();
}
