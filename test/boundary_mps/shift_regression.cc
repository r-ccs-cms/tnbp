#include "tnbp/ptns/function/bmps.h"
#include "../common/graph_fixture.h"
template<class S>
void run(const std::string& name,Edges edges,const std::vector<int>& line,
         int ownership_shift,bool zero,MPI_Comm comm) {
    run_shift_fixture<S>(name,edges,line,ownership_shift,zero,comm,
        [](const auto& l,auto& e){tnbp::bmps_shift_bonds_to_neighbors_for_edges(l,e);},
        [](auto& ctx,const auto& l,auto& e,auto& ts,auto& ids,auto& owners,MPI_Comm c){
            tnbp::bmps_shift_bonds_to_neighbors(ctx,l,e,ts,ids,owners,c);});
}
void invalid_checks(MPI_Comm comm) {
    using T=TensorOf<double>;tcapi::context_handle_t<T> ctx;tcapi::create_context(ctx);
    int rank;MPI_Comm_rank(comm,&rank);
    for(int kind=0;kind<4;++kind) {
        Edges edges{{0,2}};std::vector<int> line{0,1,2},ids;
        std::map<int,int> owners{{0,0},{1,0},{2,0}};
        std::vector<T> ts;
        if(rank==0) {
            ids=line;ts.push_back(tcapi::zeros<T>(ctx,{2}));
            ts.push_back(tcapi::zeros<T>(ctx,{}));ts.push_back(tcapi::zeros<T>(ctx,{2}));
        }
        if(kind==0) line.push_back(0);
        if(kind==1) edges.push_back({2,0});
        if(kind==2 && rank==0) {ids.pop_back();ts.pop_back();}
        if(kind==3 && rank==0) ts.back()=tcapi::zeros<T>(ctx,{3});
        const auto original=edges;const auto original_ids=ids;
        bool threw=false;
        try {tnbp::bmps_shift_bonds_to_neighbors(ctx,line,edges,ts,ids,owners,comm);}
        catch(const std::invalid_argument&) {threw=true;}
        if(!threw || edges!=original || ids!=original_ids)
            throw std::runtime_error("invalid input was not rejected before mutation");
    }
    tcapi::destroy_context(ctx);
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    try {
        invalid_checks(MPI_COMM_WORLD);
        auto cases=[&](auto scalar) {
            using S=decltype(scalar);
            for(int shift:{0,1}) {
                run<S>("triangle",{{0,1},{1,2},{0,2},{0,10},{1,11},{2,12}}, {0,1,2},shift,false,MPI_COMM_WORLD);
                run<S>("missing neighbor",{{0,2},{0,10},{1,11},{2,12}}, {0,1,2},shift,false,MPI_COMM_WORLD);
                run<S>("multiple shared neighbors",{{8,3},{8,6},{8,1},{3,6},{3,1},{6,1},{8,20},{3,21},{6,22},{1,23}},
                       {8,3,6,1},shift,false,MPI_COMM_WORLD);
                run<S>("missing neighbor with shared bonds",{{0,1},{0,2},{1,2},{0,3},{1,3},{0,20},{1,21},{2,22},{3,23}},
                       {0,1,2,3},shift,false,MPI_COMM_WORLD);
                run<S>("multi-hop",{{9,2},{9,20},{4,21},{7,22},{2,23}}, {9,4,7,2},shift,false,MPI_COMM_WORLD);
                run<S>("scalar remainder",{{0,2}}, {0,1,2},shift,false,MPI_COMM_WORLD);
                run<S>("zero",{{0,2},{1,2},{0,10},{1,11},{2,12}}, {0,1,2},shift,true,MPI_COMM_WORLD);
                run<S>("already chain",{{4,2},{2,7},{4,10},{2,11},{7,12}}, {4,2,7},shift,false,MPI_COMM_WORLD);
                Edges reversed{{0,1},{1,2},{0,2},{0,10},{1,11},{2,12}};
                std::reverse(reversed.begin(),reversed.end());
                for(auto& e:reversed) std::swap(e.first,e.second);
                run<S>("reversed edges",reversed,{0,1,2},shift,false,MPI_COMM_WORLD);
                run<S>("single site",{{0,10}}, {0},shift,false,MPI_COMM_WORLD);
                run<S>("empty line",{}, {},shift,false,MPI_COMM_WORLD);
            }
        };
        cases(double{});cases(std::complex<double>{});
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,1);
    }
    MPI_Finalize();
}
