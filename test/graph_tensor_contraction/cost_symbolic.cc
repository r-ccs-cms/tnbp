#include "tnbp/framework/graph_symbolic_region_plan.h"
#include <iostream>
#include <stdexcept>
using namespace tnbp;
void require(bool b,const char* m) {if(!b)throw std::runtime_error(m);}
int main() {
    // Nontrivial edge order; every contraction has common-neighbor fusion.
    const std::vector<std::pair<int,int>> original{{2,3},{0,2},{1,3},{0,1},{1,2},{0,3}};
    graph_string_shapes<int> strings;
    for(int site=0;site<4;++site) for(auto edge:GetSurroundingBondIndex(site,original))
        strings[site].push_back(edge%2 ? "Dx":"D");
    std::size_t cases=0;
    for(int a=0;a<4;++a)for(int b=0;b<4;++b)if(a!=b)
    for(int c=0;c<4;++c)if(c!=a && c!=b) {
        const int last=6-a-b-c;
        const std::vector<std::array<int,3>> steps{{a,b,a},{a,c,a},{a,last,a}};
        const auto symbolic=graph_schedule_cost_symbolic(original,strings,steps);
        for(unsigned d=1;d<=4;++d)for(unsigned chi=1;chi<=4;++chi) {
            auto edges=original;std::map<int,std::vector<unsigned long long>> shapes;
            for(auto entry:strings)for(auto s:entry.second)
                shapes[entry.first].push_back(detail::symbolic_power_value(detail::parse_symbolic_axis(s),d,chi));
            unsigned long long total=0;
            for(std::size_t j=0;j<steps.size();++j) {
                auto st=steps[j];const auto mac=graph_tensor_contraction_cost(edges,shapes,st[0],st[1],st[2]);total+=mac;
                require(mac==detail::symbolic_power_value(std::get<1>(symbolic)[j],d,chi),"step MAC");
                require(detail::graph_cost_elements(shapes.at(st[2]))==detail::symbolic_power_value(std::get<2>(symbolic)[j],d,chi),"step size");
                const auto prefix=graph_schedule_cost_symbolic(original,strings,
                    std::vector<std::array<int,3>>(steps.begin(),steps.begin()+j+1));
                require(edges==std::get<3>(prefix),"prefix edges");
                for(const auto& item:shapes) {
                    const auto& powers=std::get<4>(prefix).at(item.first);
                    require(item.second.size()==powers.size(),"neighbor shape rank");
                    for(std::size_t k=0;k<powers.size();++k)
                        require(item.second[k]==detail::symbolic_power_value(powers[k],d,chi),"neighbor axis dimensions");
                }
            }
            require(edges==std::get<3>(symbolic),"edge history");
            require(total==graph_symbolic_cost_value(std::get<0>(symbolic),d,chi),"total MAC");
            ++cases;
        }
    }
    auto input=strings;input.at(0)[0]="invalid";
    try {graph_contraction_cost_symbolic(original,input,0,1,0);throw std::logic_error("accepted bad syntax");}
    catch(const std::invalid_argument&) {}
    input=strings;input.at(0)[0]="1";
    try {graph_contraction_cost_symbolic(original,input,0,1,0);throw std::logic_error("accepted mismatched dimension");}
    catch(const std::invalid_argument&) {}
    auto empty=graph_schedule_cost_symbolic(original,strings,std::vector<std::array<int,3>>{});
    require(std::get<0>(empty).empty() && std::get<3>(empty)==original,"empty schedule");
    // A later invalid step must not mutate any caller input.
    const auto saved=strings;
    bool rejected=false;
    try {graph_schedule_cost_symbolic(original,strings,std::vector<std::array<int,3>>{{0,1,0},{0,1,0}});}
    catch(const std::invalid_argument&) {rejected=true;}
    require(rejected && strings==saved,"invalid schedule/input mutation");
    require(graph_symbolic_cost_value({{{2,4},1}},4,2)<graph_symbolic_cost_value({{{4,2},1}},4,2),"D>chi comparison");
    const std::map<int,int> owners{{0,2},{1,0},{2,1},{3,3}};
    auto prepared=build_graph_symbolic_region_plans(original,strings,std::vector<int>{0,1,2},owners,0,2,
        std::vector<std::array<int,3>>{{0,1,0},{0,2,0}});
    require(prepared.first.size()==3,"region tree count");
    for(const auto& plan:prepared.first) {
        require(std::get<3>(plan)==std::vector<int>({0,1,2,3}),"outside common neighbor owner missing");
        for(auto st:std::get<0>(plan))require(st[1]!=0,"target owner lost");
    }
    for(unsigned d=1;d<=4;++d)for(unsigned chi=1;chi<=4;++chi) {
        auto chosen=select_graph_symbolic_region_plan(prepared.first,prepared.second,d,chi);
        const auto& best=chosen<0 ? prepared.second:prepared.first.at(chosen);
        const auto cost=graph_symbolic_cost_value(std::get<1>(best),d,chi);
        for(const auto& plan:prepared.first)
            require(cost<=graph_symbolic_cost_value(std::get<1>(plan),d,chi),"model plan selection");
    }
    // All search branches start from the same validated input, not a previous
    // branch's fused shapes. Compare every weighted plan with standalone calls.
    for(unsigned repetitions:{1u,2u,5u}) {
        const auto plans=build_graph_symbolic_region_plans(original,strings,std::vector<int>{0,1,2,3},owners,0,repetitions,
            std::vector<std::array<int,3>>{{0,1,0},{0,2,0},{0,3,0}});
        require(plans.first.size()==18,"four-leaf enumeration changed");
        for(const auto& candidate:plans.first) {
            const auto reference=evaluate_graph_symbolic_region(original,strings,std::get<0>(candidate),owners,0,repetitions);
            require(std::get<1>(candidate)==std::get<1>(reference) && std::get<2>(candidate)==std::get<2>(reference) &&
                std::get<3>(candidate)==std::get<3>(reference) && std::get<4>(candidate)==std::get<4>(reference),"prepared search branch differs from standalone evaluation");
            for(std::size_t k=0;k<std::get<5>(candidate).size();++k)
                require(std::get<5>(candidate)[k].edges==std::get<5>(reference)[k].edges,"prepared branch graph history");
        }
        require(strings==saved,"search mutated shared shapes");
    }
    auto bad=strings;bad[99]={"D"};const auto bad_saved=bad;
    rejected=false;
    try {build_graph_symbolic_region_plans(original,bad,std::vector<int>{0,1,2},owners,0,2,
        std::vector<std::array<int,3>>{{0,1,0},{0,2,0}});}
    catch(const std::invalid_argument&){rejected=true;}
    require(rejected && bad==bad_saved,"search skipped initial shape validation or mutated input");
    require(select_graph_symbolic_region_plan(prepared.first,prepared.second,0,2)==-1,"unknown model fallback");
    std::cout<<"symbolic cost PASS: "<<cases<<" numeric schedule comparisons\n";
}
