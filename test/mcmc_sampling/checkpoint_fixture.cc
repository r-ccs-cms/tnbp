#include "tnbp/tnbp.h"
#include "../../apps/mcmc_sampling/typedef.h"
#include <filesystem>
#include <fstream>
#include <stdexcept>

// A deterministic |0,0,2> state with unsorted edges, unequal virtual dimensions,
// heterogeneous physical dimensions, and (for two ranks) an empty rank zero.
int main(int argc,char** argv) {
  if(argc!=3) return 1;
  const std::filesystem::path directory(argv[1]);
  const int ranks=std::stoi(argv[2]);
  std::filesystem::create_directories(directory);
  std::ofstream(directory/"graph.edges") << "0 2\n0 1\n";
  std::ofstream(directory/"graph.lines") << "1\n0\n2\n";
  ContextHandle ctx;tcapi::create_context(ctx);
  using Shape=tcapi::tensor_traits<Tensor>::shape_t;
  using Coors=tcapi::tensor_traits<Tensor>::elem_coors_t;
  std::vector<Tensor> tensors;
  const std::vector<Shape> shapes{{2,3,1},{3,1},{2,3}};
  const std::vector<int> states{0,0,2};
  for(std::size_t i=0;i<shapes.size();++i) {
    std::vector<Elem> data(shapes[i].back(),Elem(0));data[states[i]]=Elem(1);
    auto first=data.begin();
    tensors.push_back(tcapi::assign_from_range<Tensor>(ctx,shapes[i],first,
        [](const Coors& c) {return c.back();}));
  }
  const std::map<int,int> owners{{0,ranks-1},{1,ranks-1},{2,ranks-1}};
  for(int rank=0;rank<ranks;++rank) {
    std::ofstream out(tnbp::make_rank_filename((directory/"state").string(),rank),std::ios::binary);
    const std::vector<Tensor> empty;
    const std::vector<int> ids=rank==ranks-1 ? std::vector<int>{0,1,2} : std::vector<int>{};
    tnbp::SaveTPS(ctx,out,rank==ranks-1 ? tensors : empty,ids,owners,empty,std::vector<int>{});
    if(!out) throw std::runtime_error("checkpoint fixture write failed");
  }
  tcapi::destroy_context(ctx);
}
