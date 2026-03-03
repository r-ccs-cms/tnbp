// main.cpp
#include <iostream>
#include "tnbp/framework/graph.h"

int main(int argc, char * argv[]){

  std::vector<std::pair<int,int>> edges = {
    {0,1}, {1,2}, {2,3}, {0,4}, {1,5}, {2,6}, {3,7},
    {4,5}, {5,6}, {6,7}, {4,8}, {5,9}, {6,10}, {7,11},
    {9,10}, {10,11}, {11,12}
  };

  int i1;
  int i2;
  std::cerr << " select two integer between 0 and 12 " << std::endl;
  std::cin >> i1 >> i2;
  auto path = tnbp::FindShortestPath(edges,i1,i2);
  for(auto const & p : path) {
    std::cout << " " << p;
  }
  std::cout << std::endl;

  return 0;

}
