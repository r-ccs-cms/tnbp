#include "tnbp/tnbp.h"
#include "../common/check.h"
#include <sstream>
#include <type_traits>
template<class E> void integration(MPI_Comm comm) {
    using T=gqten::tensor<E>;using R=tcapi::real_t<T>;
    tcapi::gqten_handle ctx;tcapi::create_context(ctx);
    const std::vector<std::pair<int,int>> edges{{0,1}};
    std::vector<T> V,Emsg;std::vector<int> sites,edgeids;std::map<int,int> owners;
    tnbp::InitTensorProductState(ctx,edges,{2,2},V,sites,owners,Emsg,edgeids,comm);
    // A diagonal messenger with two known Schmidt weights exercises both
    // local and cross-rank paths of the retained-spectrum transfer.
    for(auto& tensor:V)tensor=tcapi::fill<T>(ctx,{2,2},E(1));
    for(auto& tensor:Emsg) {
        tensor=tcapi::zeros<T>(ctx,{2,2});
        tcapi::set_elem(ctx,tensor,{0,0},E(1));tcapi::set_elem(ctx,tensor,{1,1},E(0.5));
    }
    tnbp::TensorProductState<T> state(ctx,V,edges,sites,owners,comm);
    static_assert(!std::is_copy_constructible_v<tnbp::TensorProductState<T>>);
    static_assert(std::is_move_constructible_v<tnbp::TensorProductState<T>>);
    auto state_copy=state.copy(ctx);check(state_copy.NumV()==V.size(),"explicit TPS copy");
    std::vector<tcapi::bond_dim_t<T>> dims;std::vector<R> errors;
    tnbp::Truncation(ctx,edges,V,sites,owners,Emsg,edgeids,comm,1,R(0),R(0),dims,errors);
    for(auto d:dims)check(d==1,"TNBP retained dimension");
    for(auto e:errors)check(std::abs(e-R(0.2))<R(1e-5),"TNBP truncation oracle");
    for(const auto& tensor:Emsg)check(tensor.Shape()==tcapi::shape_t<T>({1,1}),"TNBP messenger remains matrix");
    std::stringstream dump(std::ios::in|std::ios::out|std::ios::binary);
    tnbp::SaveTPS(ctx,dump,V,sites,owners,Emsg,edgeids);dump.write("END",3);
    std::vector<T> loaded_v,loaded_e;std::vector<int> loaded_sites,loaded_edges;std::map<int,int> loaded_owners;
    tnbp::LoadTPS(ctx,dump,loaded_v,loaded_sites,loaded_owners,loaded_e,loaded_edges);
    check(sites==loaded_sites && edgeids==loaded_edges && owners==loaded_owners,"TPS stream metadata");
    for(std::size_t i=0;i<V.size();++i)check(tcapi::close(ctx,V[i],loaded_v[i],R(0)),"TPS stream vertex");
    for(std::size_t i=0;i<Emsg.size();++i)check(tcapi::close(ctx,Emsg[i],loaded_e[i],R(0)),"TPS stream messenger");
    char suffix[3];dump.read(suffix,3);check(std::string(suffix,3)=="END","TPS stream trailing bytes");
    auto matrix=tcapi::zeros<T>(ctx,{2,2});tcapi::set_elem(ctx,matrix,{0,0},E(4));tcapi::set_elem(ctx,matrix,{1,1},E(9));T root;
    tnbp::SquareRoot(ctx,matrix,1,root);
    check(std::abs(root.GetElem(0)-E(2))<R(1e-5) && std::abs(root.GetElem(3)-E(3))<R(1e-5),"SquareRoot migrated explicit labels");
    tcapi::destroy_context(ctx);
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    try {
        integration<double>(MPI_COMM_WORLD);integration<std::complex<double>>(MPI_COMM_WORLD);
        std::cout<<"PASS "<<checks<<" checks per rank\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
    MPI_Finalize();return 0;
}
