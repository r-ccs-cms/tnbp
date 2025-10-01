/// This file is a part of r-ccs-cms/tnbp
/**
 */
#ifndef TNBP_FRAMEWORK_GRAPH_H
#define TNBP_FRAMEWORK_GRAPH_H

#include <vector>
#include <queue>

#include <vector>
#include <utility>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
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

  template <typename IntT>
  std::vector<std::vector<std::pair<IntT,IntT>>> GreedyEdgeLayering(
			     std::vector<std::pair<IntT,IntT>> edges) {
    
    std::vector<std::vector<std::pair<IntT,IntT>>> layers;
    
    // set of temporary edges
    std::unordered_set<size_t> used_indices;
    std::vector<bool> used(edges.size(), false);
    
    while (true) {
      std::unordered_set<IntT> used_vertices;
      std::vector<std::pair<IntT,IntT>> layer;
      
      for (size_t i = 0; i < edges.size(); ++i) {
	if (used[i]) continue;
	
	auto [u, v] = edges[i];
	if (used_vertices.count(u) || used_vertices.count(v)) continue;
	
	// if it is not used, add it to layer
	layer.push_back(edges[i]);
	used_vertices.insert(u);
	used_vertices.insert(v);
	used[i] = true;
      }
      
      if (layer.empty()) break; // if layer was empy quit
      layers.push_back(std::move(layer));
    }
    
    return layers;
  }

  template <std::integral IntT>
  std::vector<std::vector<IntT>> BuildContractionLayers(
	 const std::vector<std::pair<IntT, IntT>>& edges,
	 const std::vector<IntT>& initial_boundary
							) {
    
    using NeighborMap = std::unordered_map<IntT, std::vector<IntT>>;

    NeighborMap neighbors;
    for (const auto& [u, v] : edges) {
      neighbors[u].push_back(v);
      neighbors[v].push_back(u);
    }
    
    std::unordered_set<IntT> absorbed(initial_boundary.begin(), initial_boundary.end());
    
    std::vector<std::vector<IntT>> layers;
    layers.push_back(initial_boundary);  //
    
    std::set<IntT> all_nodes;
    for (const auto& [u, v] : edges) {
      all_nodes.insert(u);
      all_nodes.insert(v);
    }
    
    while (true) {
      std::vector<IntT> next_layer;
      for (const auto& node : all_nodes) {
	if (absorbed.count(node)) continue;
	
	for (const auto& neighbor : neighbors[node]) {
	  if (absorbed.count(neighbor)) {
	    next_layer.push_back(node);
	    break;
	  }
	}
      }
      
      if (next_layer.empty()) break;
      
      std::sort(next_layer.begin(), next_layer.end());
      next_layer.erase(std::unique(next_layer.begin(), next_layer.end()), next_layer.end());
      
      for (const auto& node : next_layer) {
	absorbed.insert(node);
      }
      layers.push_back(std::move(next_layer));
    }
    
    return layers;
  }


/**
 * @brief Extract edges from a global edge list that are fully contained in a given line.
 *
 * This function takes a set of edges representing the global graph and a list of vertices
 * representing a single line. It returns only those edges whose both endpoints are included
 * in the specified line.
 *
 * @tparam IntT Integral type for vertex labels.
 * @param global_edge The full set of edges in the graph, represented as pairs of vertex labels.
 * @param line A list of vertex labels defining the line.
 * @return std::vector<std::pair<IntT, IntT>> The subset of edges contained entirely within the line.
 *
 * @note Vertex labels are kept as in the original global graph (no reindexing).
 *
 * @complexity O(E) where E is the number of edges in global_edge.
 */
  template <std::integral IntT>
  std::vector<std::pair<IntT, IntT>>
  ExtractLineEdges(const std::vector<std::pair<IntT, IntT>>& global_edge,
		   const std::vector<IntT>& line) {
    std::unordered_set<IntT> line_set(line.begin(), line.end());
    
    std::vector<std::pair<IntT, IntT>> line_edge;
    line_edge.reserve(global_edge.size());
    
    for (auto& [u, v] : global_edge) {
      if (line_set.count(u) && line_set.count(v)) {
	line_edge.emplace_back(u, v);
      }
    }
    return line_edge;
  }
  
}

#endif

