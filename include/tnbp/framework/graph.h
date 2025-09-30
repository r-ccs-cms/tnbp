/// This file is a part of r-ccs-cms/tnbp
/**
 */
#ifndef TNBP_FRAMEWORK_GRAPH_H
#define TNBP_FRAMEWORK_GRAPH_H

#include <vector>
#include <queue>

#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <limits>
#include <concepts>

namespace tnbp {

  /**
     Function to extract sites from edges/bonds
   */
  template <std::integral IntT>
  std::vector<IntT> GetSiteIndexFromBond(
	const std::vector<std::pair<IntT,IntT>> edge) {
    std::vector<IntT> site;
    site.reserve(edge.size()*2);
    for(auto [u,v] : edge) {
      site.push_back(u);
      site.push_back(v);
    }
    std::sort(site.begin(),site.end());
    site.erase(std::unique(site.begin(),site.end()),site.end());
    return site;
  }

  /**
     Function to extract surrounding edge labels
   */
  template <std::integral IntT>
  std::vector<IntT> GetSurroundingBondIndex(IntT SiteIdx,
					    const std::vector<std::pair<IntT,IntT>> & edge) {
    std::vector<IntT> res;
    IntT m=0;
    for(auto const & [i,j] : edge) {
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
  
  template <std::integral IntT>
  std::vector<IntT> FindShortestPath(const std::vector<std::pair<IntT,IntT>> & I,
                                    IntT i1, IntT i2) {
    std::vector<IntT> Q = GetSiteIndexFromBond(I);
    std::vector<std::vector<IntT>> adj(Q.size());
    for(const auto & edge : I) {
      adj[edge.first].push_back(edge.second);
      adj[edge.second].push_back(edge.first);
    }
    std::queue<IntT> q;
    std::vector<bool> visited(Q.size(),false);
    std::vector<IntT> prev(Q.size(),-1);

    q.push(i1);
    visited[i1] = true;

    while ( !q.empty() ) {
      IntT node = q.front();
      q.pop();
      if( node == i2 ) {
        std::vector<IntT> path;
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

  template <std::integral IntT>
  std::vector<std::pair<IntT,IntT>> CompressVertices(
	       std::vector<std::pair<IntT,IntT>> edge) {
    std::vector<IntT> site;
    site.reserve(edge.size()*2);
    for(auto [u,v] : edge) {
      site.push_back(u);
      site.push_back(v);
    }
    std::sort(site.begin(),site.end());
    site.erase(std::unique(site.begin(),site.end()),site.end());

    if(site.size() > static_cast<std::size_t>(std::numeric_limits<IntT>::max())) {
      throw std::overflow_error("CompressVertices: new label exceeds IntT range");
    }

    std::unordered_map<IntT,IntT> to_new;
    for(std::size_t i=0; i < site.size(); i++) {
      to_new.try_emplace(site[i],static_cast<IntT>(i));
    }

    for(auto & [u,v] : edge) {
      u = to_new[u];
      v = to_new[v];
    }
    return edge;
  }
  
  
}

#endif

