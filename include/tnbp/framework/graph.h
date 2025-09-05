/// This file is a part of r-ccs-cms/tnbp
/**
 */
#ifndef TNBP_FRAMEWORK_GRAPH_H
#define TNBP_FRAMEWORK_GRAPH_H

#include <vector>
#include <queue>

namespace tnbp {

  std::vector<int> GetSiteIndexFromBond(
		      const std::vector<std::pair<int,int>> & I) {
    std::vector<int> res(2*I.size());
    auto itres = res.begin();
    for(auto const & [i,j] : I) {
      *itres++ = i;
      *itres++ = j;
    }
    std::sort(res.begin(),res.end());
    res.erase(std::unique(res.begin(),res.end()),res.end());
    return res;
  }

  std::vector<int> GetSurroundingBondIndex(int SiteIdx,
		      const std::vector<std::pair<int,int>> & I) {
    std::vector<int> res;
    int m=0;
    for(auto const & [i,j] : I) {
      if( i == SiteIdx ) {
	res.push_back(m);
      }
      if( j == SiteIdx ) {
	res.push_back(m);
      }
      m++;
    }
    return res;
  }

  std::vector<int> FindShortestPath(const std::vector<std::pair<int,int>> & I,
                                    int i1, int i2) {
    std::vector<int> Q = GetSiteIndexFromBond(I);
    std::vector<std::vector<int>> adj(Q.size());
    for(const auto & edge : I) {
      adj[edge.first].push_back(edge.second);
      adj[edge.second].push_back(edge.first);
    }
    std::queue<int> q;
    std::vector<bool> visited(Q.size(),false);
    std::vector<int> prev(Q.size(),-1);

    q.push(i1);
    visited[i1] = true;

    while ( !q.empty() ) {
      int node = q.front();
      q.pop();
      if( node == i2 ) {
        std::vector<int> path;
        for(int at = i2; at != -1; at =prev[at]) {
          path.push_back(at);
        }
        std::reverse(path.begin(),path.end());
        return path;
      }

      for(int neighbor : adj[node]) {
        if (!visited[neighbor]) {
          q.push(neighbor);
          visited[neighbor] = true;
          prev[neighbor] = node;
        }
      }
    }
    return {};
  }
  
  
}

#endif

