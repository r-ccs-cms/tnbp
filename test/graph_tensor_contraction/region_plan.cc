#include "tnbp/framework/graph_region_contraction_plan.h"
#include <iostream>
using Axes=std::vector<std::vector<int>>;
using Dims=std::vector<std::uint64_t>;
using Plan=std::vector<std::array<std::size_t,3>>;
void require(bool b,const char* m){if(!b)throw std::runtime_error(m);}
// Independent exhaustive oracle: leaf subsets and original edges, with no
// planner metadata, graph update helpers or binary-join evaluator.
std::pair<std::uint64_t,std::uint64_t> oracle(const std::vector<std::pair<int,int>>& edges,
    const Dims& dims,int n,int p) {
    std::pair<std::uint64_t,std::uint64_t> best{UINT64_MAX,UINT64_MAX};
    std::vector<unsigned> initial;for(int i=0;i<n;++i) initial.push_back(1u<<i);
    auto search=[&](auto&& self,std::vector<unsigned> live,std::uint64_t total,std::uint64_t peak)->void {
        if(live.size()==1) {best=std::min(best,std::make_pair(total,peak));return;}
        for(std::size_t i=0;i<live.size();++i)for(std::size_t j=i+1;j<live.size();++j){
            unsigned a=live[i],b=live[j],c=a|b;bool connected=false;std::uint64_t cost=1,size=1;
            for(std::size_t k=0;k<edges.size();++k){
                const auto [x,y]=edges[k];
                const bool ax=bool(a&(1u<<x)),ay=bool(a&(1u<<y));
                const bool bx=bool(b&(1u<<x)),by=bool(b&(1u<<y));
                if((ax!=ay)||(bx!=by))cost*=dims[k];
                if((ax&&by)||(ay&&bx))connected=true;
                if(bool(c&(1u<<x))!=bool(c&(1u<<y)))size*=dims[k];
            }
            if(!connected)continue;
            auto next=live;next.erase(next.begin()+j);next.erase(next.begin()+i);next.push_back(c);
            self(self,next,total+cost*((c&1)?p:1),std::max(peak,size));
        }
    };search(search,initial,0,0);return best;
}
Plan baseline_chain(int n){Plan s;std::size_t current=0;for(int i=1;i<n;++i){s.push_back({current,std::size_t(i),std::size_t(n+i-1)});current=n+i-1;}return s;}
int main(){try{
    int cases=0;
    for(int n=1;n<=5;++n)for(int complete:{0,1})for(int p:{1,2,5}) {
        std::vector<std::pair<int,int>> edges;Dims dims;Axes axes(n);
        for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if(complete||j==i+1){
            const int k=edges.size();edges.emplace_back(i,j);dims.push_back(2+(k%2));axes[i].push_back(k);axes[j].push_back(k);
        }
        Plan plan;std::uint64_t peak;bool fallback;
        auto cost=tnbp::make_small_region_contraction_plan(axes,dims,std::vector<int>(n,0),0,p,baseline_chain(n),plan,peak,fallback);
        require(!fallback && std::make_pair(cost,peak)==oracle(edges,dims,n,p),"exhaustive subset oracle mismatch");
        Plan again;std::uint64_t peak2;bool fallback2;
        auto cost2=tnbp::make_small_region_contraction_plan(axes,dims,std::vector<int>(n,0),0,p,baseline_chain(n),again,peak2,fallback2);
        require(plan==again && cost==cost2 && peak==peak2,"non-deterministic plan");++cases;
    }
    // Triangle: baseline total overflows; eliminating the large shared edge
    // first keeps the total representable (large+1).
    const auto large=UINT64_MAX/2+1;
    Axes triangle{{0,1},{0,2},{1,2}};Plan plan;std::uint64_t peak=77;bool fallback=true;
    auto cost=tnbp::make_small_region_contraction_plan(triangle,Dims{large,1,1},std::vector<int>{0,1,2},0,1,
        Plan{{1,2,3},{0,3,4}},plan,peak,fallback);
    require(cost==large+1 && !fallback,"triangle optimal cost");
    // Candidate weighting can favor constructing the whole common subtree.
    Axes square{{0,3,4,6},{0,1,2},{1,3,5},{2,4,7}};Dims dims{2,3,3,2,2,3,2,3};
    Plan old{{2,1,4},{0,4,5},{5,3,6}};
    cost=tnbp::make_small_region_contraction_plan(square,dims,std::vector<int>{2,0,1,1},0,2,old,plan,peak,fallback);
    require(cost==612 && peak==36 && plan[0][0]==1,"square MAC/owner normalization");
    cost=tnbp::make_small_region_contraction_plan(square,dims,std::vector<int>{2,0,1,1},0,3,old,plan,peak,fallback);
    require(cost==756 && peak==72,"candidate common-subtree weighting");
    // Six-site fallback keeps legacy orientation, even when min leaf differs.
    Axes chain(6);Dims cd(5,2);for(int i=0;i<5;++i){chain[i].push_back(i);chain[i+1].push_back(i);}
    Plan old6{{2,1,6},{0,6,7},{7,3,8},{8,4,9},{9,5,10}};
    tnbp::make_small_region_contraction_plan(chain,cd,std::vector<int>{0,1,2,3,4,5},0,2,old6,plan,peak,fallback);
    require(fallback && plan.front()[0]==2,"fallback owner orientation changed");
    for(int kind=0;kind<4;++kind){
        auto bad_axes=triangle;auto bad_dims=Dims{2,3,5};auto oldplan=plan;auto oldpeak=peak;auto oldfallback=fallback;bool caught=false;
        if(kind==0)bad_axes[0].push_back(0);
        if(kind==1)bad_dims[0]=0;
        if(kind==2)bad_dims={UINT64_MAX,2,2};
        if(kind==3)bad_dims={2,3,5};
        try{tnbp::make_small_region_contraction_plan(bad_axes,bad_dims,std::vector<int>{0,1,2},0,
            kind==3?UINT64_MAX:2,baseline_chain(3),plan,peak,fallback);}catch(const std::exception&){caught=true;}
        require(caught && plan==oldplan && peak==oldpeak && fallback==oldfallback,"planner failure not atomic");
    }
    std::cout<<cases<<" independent exhaustive oracle cases; reuse, fallback, owner, overflow and atomicity PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
