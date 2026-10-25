#include "tnbp/framework/graph_contraction_cost.h"
#include <iostream>
#include <string>
using Edges=std::vector<std::pair<int,int>>;
using Shapes=std::map<int,std::vector<std::uint64_t>>;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
Shapes shapes_for(const Edges& edges,const std::vector<std::uint64_t>& dims){
    Shapes shapes;
    for(auto site:tnbp::GetSiteIndexFromBond(edges))
        for(auto e:tnbp::GetSurroundingBondIndex(site,edges)) shapes[site].push_back(dims.at(e));
    return shapes;
}
template<class F> void reject(F call){bool caught=false;try{call();}catch(const std::invalid_argument&){caught=true;}require(caught,"invalid input accepted");}
void analytic_checks(){
    {
        Edges e{{0,1},{1,2}};auto s=shapes_for(e,{2,3});const auto schedule=e;
        std::vector<std::uint64_t> costs,sizes;std::vector<std::array<int,3>> steps;
        auto total=tnbp::graph_tensor_network_contraction_cost(e,s,schedule,{}, {},&costs,&sizes,&steps);
        require(total==9 && costs==std::vector<std::uint64_t>({6,3}) && sizes==std::vector<std::uint64_t>({3,1}),"chain arithmetic");
        require(e.empty() && s==Shapes({{0,{}}}) && steps==std::vector<std::array<int,3>>({{0,1,0},{0,2,0}}),"chain final metadata");
    }
    const Edges triangle{{0,1},{0,2},{1,2}};
    for(auto first:triangle){
        auto e=triangle;auto s=shapes_for(e,{2,3,5});Edges schedule{first};
        for(auto edge:triangle)if(edge!=first)schedule.push_back(edge);
        const auto total=tnbp::graph_tensor_network_contraction_cost(e,s,schedule);
        const auto expected=first==triangle[0]?45u:first==triangle[1]?40u:36u;
        require(total==expected && e.empty() && s.size()==1,"triangle cost depends on order");
    }
    {
        Edges e{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};auto s=shapes_for(e,{2,3,5,7,11,13});auto schedule=e;
        std::vector<std::uint64_t> costs,sizes;std::vector<Shapes> history;
        auto total=tnbp::graph_tensor_network_contraction_cost(e,s,schedule,{}, {},&costs,&sizes,nullptr,&history);
        require(total==18040 && costs==std::vector<std::uint64_t>({2310,15015,715}),"K4 MACs");
        require(sizes==std::vector<std::uint64_t>({1155,715,1}),"K4 result sizes");
        require(history[0]==Shapes({{0,{21,55}},{2,{21,13}},{3,{55,13}}}),"K4 updates both outside neighbors");
    }
    {
        Edges e{{0,1},{0,2},{0,3},{1,2},{2,3},{1,4},{2,5},{3,6}};
        auto s=shapes_for(e,{3,2,3,2,2,3,2,3});std::vector<std::uint64_t> costs,sizes;std::vector<Shapes> history;
        const auto total=tnbp::graph_tensor_network_contraction_cost(e,s,Edges{{1,0},{2,1},{2,3}},{{1,1},{2,2}},{},&costs,&sizes,nullptr,&history);
        require(total==360 && costs==std::vector<std::uint64_t>({108,144,108}) && sizes==std::vector<std::uint64_t>({36,36,18}),"square D2 chi3 arithmetic");
        require(e==Edges({{2,4},{2,5},{2,6}}) && s.at(2)==std::vector<std::uint64_t>({3,2,3}),"square surviving edge/axis order");
        require(history[0].at(2)==std::vector<std::uint64_t>({4,2,2}) && history[1].at(3)==std::vector<std::uint64_t>({6,3}),"square target and boundary fusion");
    }
    {
        auto e=triangle;auto s=shapes_for(e,{2,3,5});std::vector<std::array<int,3>> steps;
        auto total=tnbp::graph_tensor_network_contraction_cost(e,s,triangle,{{1,5}},{{1,-9}},nullptr,nullptr,&steps);
        require(total==45 && s==Shapes({{-9,{}}}) && steps==std::vector<std::array<int,3>>({{1,0,-9},{-9,2,-9}}),"priority/label inheritance");
        e=triangle;s=shapes_for(e,{2,3,5});
        total=tnbp::graph_tensor_network_contraction_cost(e,s,Edges{{2,1},{0,1},{0,2}}, {}, {},nullptr,nullptr,&steps);
        require(total==36 && steps[0]==std::array<int,3>({2,1,2}),"tie follows scheduled endpoint orientation");
    }
    {
        auto e=triangle;auto s=shapes_for(e,{2,3,5});s[99]={};const auto original=s;
        std::vector<std::uint64_t> costs{99};
        require(tnbp::graph_tensor_network_contraction_cost(e,s,Edges{}, {}, {},&costs)==0 && costs.empty() && s==original && e==triangle,"empty schedule preserves input and clears output");
        require(tnbp::graph_tensor_network_contraction_cost(e,s,Edges{{0,1}})==30 && e==Edges({{0,2}}) && s.at(2)==std::vector<std::uint64_t>({15}) && s.count(99),"partial schedule and isolated scalar");
    }
    std::cout<<"analytic chain/triangle/K4/square, outside shapes, priorities/labels, partial/empty PASS\n";
}
void schedule_permutations(){
    const Edges graph{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    const auto initial=shapes_for(graph,{2,3,5,7,11,13});
    std::vector<int> order{0,1,2,3,4,5};std::set<std::uint64_t> totals;int checked=0;
    auto count=[](const std::vector<std::uint64_t>& shape){std::uint64_t n=1;for(auto d:shape)n*=d;return n;};
    do{
        Edges schedule;for(auto i:order)schedule.push_back(graph[i]);
        auto edges=graph;auto shapes=initial;
        std::vector<std::uint64_t> costs;std::vector<std::array<int,3>> steps;
        std::vector<Shapes> shapes_trace;std::vector<Edges> edges_trace;
        const auto total=tnbp::graph_tensor_network_contraction_cost(edges,shapes,schedule,{}, {},&costs,nullptr,&steps,&shapes_trace,&edges_trace);
        require(edges.empty() && shapes.size()==1 && shapes.begin()->second.empty() && steps.size()==3,"full K4 schedule did not close");
        std::uint64_t reference_total=0;
        for(std::size_t i=0;i<steps.size();++i){
            const auto& before_shapes=i?shapes_trace[i-1]:initial;
            const auto& before_edges=i?edges_trace[i-1]:graph;
            const auto a=steps[i][0],b=steps[i][1];
            const auto edge=tnbp::find_edge(before_edges,a,b)-before_edges.begin();
            const auto bonds=tnbp::GetSurroundingBondIndex(a,before_edges);
            const auto axis=std::find(bonds.begin(),bonds.end(),edge)-bonds.begin();
            // Independent input-size/shared-dimension formula, for small exact integers.
            const auto expected=count(before_shapes.at(a))*count(before_shapes.at(b))/before_shapes.at(a).at(axis);
            require(costs[i]==expected,"MAC differs from input-product/shared-dimension formula");reference_total+=expected;
        }
        require(total==reference_total,"permutation total MACs");totals.insert(total);++checked;
    }while(std::next_permutation(order.begin(),order.end()));
    require(checked==720 && totals.size()>1,"schedule permutation coverage");
    std::cout<<"720 K4 schedule permutations and independent per-step MAC formula PASS\n";
}
void invalid_checks(){
    const Edges original{{0,1},{0,2},{1,2}};const auto initial=shapes_for(original,{2,3,5});
    for(int kind=0;kind<9;++kind){
        auto e=original;auto s=initial;Edges schedule=e;
        if(kind==0)e.push_back({1,0});
        if(kind==1)e.push_back({2,2});
        if(kind==2)s.erase(2);
        if(kind==3)s[0].pop_back();
        if(kind==4)s[0][0]=0;
        if(kind==5)s[1][0]=3;
        if(kind==6)s[99]={2};
        if(kind==7)schedule={{0,99}};
        if(kind==8)schedule.push_back(schedule.front());
        const auto before_e=e;const auto before_s=s;std::vector<std::uint64_t> costs{99},sizes{88};
        reject([&]{tnbp::graph_tensor_network_contraction_cost(e,s,schedule,{}, {},&costs,&sizes);});
        require(e==before_e && s==before_s && costs==std::vector<std::uint64_t>({99}) && sizes==std::vector<std::uint64_t>({88}),"failure modified input/output");
    }
    for(int result:{1,2,99}){
        auto e=original;auto s=initial;s[99]={};
        reject([&]{tnbp::graph_tensor_contraction_cost(e,s,0,1,result);});
    }
    {
        Edges e{{0,1}};std::map<int,std::vector<int>> s{{0,{-2}},{1,{-2}}};reject([&]{tnbp::graph_tensor_contraction_cost(e,s,0,1,0);});
        auto u=shapes_for(e,{2});std::vector<std::uint64_t> out;
        reject([&]{tnbp::graph_tensor_network_contraction_cost(e,u,e,{}, {},&out,&out);});
    }
    for(int kind=0;kind<2;++kind){
        Edges e=kind==0?Edges{{0,1},{0,2},{1,3}}:Edges{{0,1},{2,3}};
        const auto big=std::numeric_limits<std::uint64_t>::max();
        auto s=shapes_for(e,kind==0?std::vector<std::uint64_t>{1,big,2}:std::vector<std::uint64_t>{big/2+1,big/2+1});
        const auto before=e;const auto old=s;std::vector<std::uint64_t> out{77};bool caught=false;
        try{tnbp::graph_tensor_network_contraction_cost(e,s,kind==0?Edges{{0,1}}:e,{}, {},&out);}catch(const std::overflow_error&){caught=true;}
        require(caught && e==before && s==old && out==std::vector<std::uint64_t>({77}),"MAC overflow not atomic");
    }
    {
        auto e=original;std::map<int,std::vector<unsigned>> s{{0,{1,65536}},{1,{1,65536}},{2,{65536,65536}}};const auto old=s;
        bool caught=false;try{tnbp::graph_tensor_contraction_cost(e,s,0,1,0);}catch(const std::overflow_error&){caught=true;}
        require(caught && e==original && s==old,"dimension type overflow");
    }
    {
        // First neighbor is already reformed when the second fusion overflows.
        Edges e{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
        std::map<int,std::vector<unsigned>> s;
        const std::vector<unsigned> dims{1,2,65536,3,65536,1};
        for(auto site:tnbp::GetSiteIndexFromBond(e))
            for(auto bond:tnbp::GetSurroundingBondIndex(site,e)) s[site].push_back(dims[bond]);
        const auto before_e=e;const auto before_s=s;bool caught=false;
        try{tnbp::graph_tensor_contraction_cost(e,s,0,1,0);}
        catch(const std::overflow_error&){caught=true;}
        require(caught && e==before_e && s==before_s,"partial neighbor mutation escaped single-step failure");
    }
    {
        // A successful step and populated private reports precede this failure.
        auto e=original;auto s=initial;std::vector<std::uint64_t> costs{99},sizes{88};
        std::vector<std::array<int,3>> steps{{9,8,7}};
        std::vector<Shapes> history{initial};std::vector<Edges> edges_history{original};
        reject([&]{tnbp::graph_tensor_network_contraction_cost(e,s,Edges{{0,1},{0,99}},
            {}, {},&costs,&sizes,&steps,&history,&edges_history);});
        require(e==original && s==initial && costs==std::vector<std::uint64_t>{99} &&
            sizes==std::vector<std::uint64_t>{88} && steps==std::vector<std::array<int,3>>{{9,8,7}} &&
            history==std::vector<Shapes>{initial} && edges_history==std::vector<Edges>{original},
            "late schedule failure modified inputs or reports");
    }
    std::cout<<"validation, overflow, exception atomicity PASS\n";
}
int main(){try{analytic_checks();schedule_permutations();invalid_checks();std::cout<<"shape-only graph cost PASS (no tensor backend/MPI)\n";return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
