/// This file is a part of r-ccs-cms/tnbp
/**
 */
#ifndef TNBP_FRAMEWORK_GRAPH_H
#define TNBP_FRAMEWORK_GRAPH_H

#include <vector>
#include <queue>
#include <utility>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <limits>

namespace tnbp {

  /**
     Definition of regulated edges
  */
  template <typename IntT>
  std::pair<IntT, IntT> make_edge(IntT u, IntT v) {
    if (u > v) std::swap(u, v);
    return {u, v};
  }
  
  template <typename IntT>
  constexpr bool same_edge(const std::pair<IntT,IntT>& a, const std::pair<IntT,IntT>& b) {
    return make_edge(a.first, a.second) == make_edge(b.first, b.second);
  }

  template <typename IntT>
  auto find_edge(std::vector<std::pair<IntT,IntT>>& edges, IntT u, IntT v) {
    std::pair<IntT,IntT> target = make_edge(u, v);
    return std::find_if(edges.begin(), edges.end(),
        [&](const auto& e) {
            return make_edge(e.first, e.second) == target;
        });
  }

  template <typename IntT>
  auto find_edge(const std::vector<std::pair<IntT,IntT>>& edges, IntT u, IntT v) {
    std::pair<IntT,IntT> target = make_edge(u, v);
    return std::find_if(edges.begin(), edges.end(),
			[&](const auto& e) {
			  return make_edge(e.first, e.second) == target;
			});
  }

  template <typename IntT>
  bool contains_edge(const std::vector<std::pair<IntT,IntT>>& edges, IntT u, IntT v) {
    return find_edge(edges, u, v) != edges.end();
  }

  template <typename IntT>
  bool add_edge(std::vector<std::pair<IntT,IntT>>& edges, IntT u, IntT v) {
    edges.push_back(make_edge(u, v));
    return true;
  }

  template <typename IntT>
  bool add_edge_if_absent(std::vector<std::pair<IntT,IntT>>& edges, IntT u, IntT v) {
    std::pair<IntT,IntT> e = make_edge(u, v);
    if (std::find(edges.begin(), edges.end(), e) != edges.end()) return false;
    edges.push_back(e);
    return true;
  }

  template <typename IntT>
  bool erase_edge(std::vector<std::pair<IntT,IntT>>& edges, IntT u, IntT v) {
    auto it = find_edge(edges, u, v);
    if (it == edges.end()) return false;
    edges.erase(it);
    return true;
  }

  template <typename IntT>
  std::size_t erase_all_edges(std::vector<std::pair<IntT,IntT>>& edges, IntT u, IntT v) {
    std::pair<IntT,IntT> e = make_edge(u, v);
    auto old_size = edges.size();
    edges.erase(std::remove(edges.begin(), edges.end(), e), edges.end());
    return old_size - edges.size();
  }

  template <typename IntT>
  constexpr void normalize_edge(std::pair<IntT, IntT>& e) {
    if (e.second < e.first) std::swap(e.first, e.second);
  }

  template <typename IntT>
  void normalize_edges(std::vector<std::pair<IntT, IntT>>& edges) {
    for (auto& e : edges) {
      normalize_edge(e);
    }
  }

  template <typename IntT>
  void normalize_and_unique_edges(std::vector<std::pair<IntT, IntT>>& edges) {
    for (auto& e : edges) {
      normalize_edge(e);
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  }
  
  
  /**
     Function to extract sites from edges/bonds
   */
  template <typename IntT>
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
  template <typename IntT>
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

  template <typename IntT>
  std::vector<IntT> FindShortestPath(const std::vector<std::pair<IntT,IntT>> &E,
				     IntT src, IntT dst)
  {
    // 1. Compress labels
    std::vector<IntT> sites;
    sites.reserve(E.size()*2);
    for(auto [u,v] : E){ sites.push_back(u); sites.push_back(v); }
    std::sort(sites.begin(), sites.end());
    sites.erase(std::unique(sites.begin(), sites.end()), sites.end());
    
    std::unordered_map<IntT,std::size_t> to_idx;
    to_idx.reserve(sites.size()*2);
    for(std::size_t i=0;i<sites.size();++i) to_idx[sites[i]]=i;
    
    if(!to_idx.count(src) || !to_idx.count(dst)) return {};
    
    std::size_t s = to_idx[src], t = to_idx[dst];
    std::vector<std::vector<std::size_t>> adj(sites.size());
    for(auto [u,v]:E){
      adj[to_idx[u]].push_back(to_idx[v]);
      adj[to_idx[v]].push_back(to_idx[u]);
    }
    
    // 2. BFS
    std::queue<std::size_t> q;
    std::vector<char> visited(sites.size(),0);
    std::vector<std::size_t> prev(sites.size(),-1);
    
    q.push(s); visited[s]=1;
    while(!q.empty()){
      auto u=q.front(); q.pop();
      if(u==t) break;
      for(auto v:adj[u]){
	if(!visited[v]){
	  visited[v]=1;
	  prev[v]=u;
	  q.push(v);
	}
      }
    }
    
    if(!visited[t]) return {};
    
    // 3. Reconstruction of path
    std::vector<IntT> path;
    for(std::size_t cur=t;cur!=(std::size_t)-1;cur=prev[cur]){
      path.push_back(sites[cur]);
    }
    std::reverse(path.begin(),path.end());
    return path;
  }

  template <typename IntT>
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

  template <typename IntT>
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
  template <typename IntT>
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

/**
 * @brief Extract edges that connect two disjoint vertex sets (line_a, line_b) from a global edge list.
 *
 * Given a global edge list and two vertex lists line_a and line_b, this function returns
 * only those edges with one endpoint in line_a and the other in line_b.
 *
 * @tparam IntT Integral type for vertex labels.
 * @param global_edges The full set of edges in the graph (undirected), as pairs of vertex labels.
 * @param line_a Vertex labels for layer A.
 * @param line_b Vertex labels for layer B.
 * @return std::vector<std::pair<IntT,IntT>> Edges that connect line_a and line_b (original labels).
 *
 * @note An edge (u,v) is included iff (u ∈ A and v ∈ B) or (u ∈ B and v ∈ A).
 * @complexity O(E) where E = global_edges.size().
 */
  template <typename IntT>
  std::vector<std::pair<IntT,IntT>>
  ExtractInterlayerEdges(const std::vector<std::pair<IntT,IntT>>& global_edges,
			 const std::vector<IntT>& line_a,
			 const std::vector<IntT>& line_b) {
    std::unordered_set<IntT> A(line_a.begin(), line_a.end());
    std::unordered_set<IntT> B(line_b.begin(), line_b.end());
    
    std::vector<std::pair<IntT,IntT>> inter_edges;
    inter_edges.reserve(global_edges.size());
    
    for (const auto& [u,v] : global_edges) {
      const bool uA = A.count(u), vA = A.count(v);
      const bool uB = B.count(u), vB = B.count(v);
      if ( (uA && vB) || (uB && vA) ) {
	inter_edges.emplace_back(u, v);
      }
    }
    return inter_edges;
  }
  
/**
 * @brief Extract interlayer edges and reindex endpoints locally per layer.
 *
 * Vertices in line_a are mapped to [0, |line_a|), and vertices in line_b are mapped to [0, |line_b|).
 * Returned edges keep the original side information (A/B) only by where endpoints came from; the pairs
 * are NOT forced into a single orientation.
 *
 * @tparam IntT Integral type for vertex labels.
 * @param global_edges The full set of edges in the graph (undirected), as pairs of vertex labels.
 * @param line_a Vertex labels for layer A (reindexed to [0, |A|)).
 * @param line_b Vertex labels for layer B (reindexed to [0, |B|)).
 * @return std::vector<std::pair<IntT,IntT>> Interlayer edges with endpoints reindexed
 *         (A-side in [0,|A|), B-side in [0,|B|)). The first/second of the pair preserves the original ordering.
 *
 * @note If you need a canonical orientation (A->B), use the oriented variant below.
 * @complexity O(E) where E = global_edges.size().
 */
  template <typename IntT>
  std::vector<std::pair<IntT,IntT>>
  ExtractInterlayerEdgesReindexed(const std::vector<std::pair<IntT,IntT>>& global_edges,
				  const std::vector<IntT>& line_a,
				  const std::vector<IntT>& line_b) {
    std::unordered_map<IntT, IntT> a_local, b_local;
    a_local.reserve(line_a.size());
    b_local.reserve(line_b.size());
    
    for (size_t i = 0; i < line_a.size(); ++i) a_local[line_a[i]] = static_cast<IntT>(i);
    for (size_t i = 0; i < line_b.size(); ++i) b_local[line_b[i]] = static_cast<IntT>(i);
    
    std::vector<std::pair<IntT,IntT>> inter_edges;
    inter_edges.reserve(global_edges.size());
    
    for (const auto& [u,v] : global_edges) {
      auto iuA = a_local.find(u), ivA = a_local.find(v);
      auto iuB = b_local.find(u), ivB = b_local.find(v);
      
      if (iuA != a_local.end() && ivB != b_local.end()) {
	// u in A, v in B
	inter_edges.emplace_back(iuA->second, ivB->second);
      } else if (iuB != b_local.end() && ivA != a_local.end()) {
	// u in B, v in A
	inter_edges.emplace_back(ivA->second, iuB->second);
      }
    }
    return inter_edges;
  }

/**
 * @brief Extract interlayer edges oriented from A to B without reindexing.
 *
 * Given an undirected global edge list and two disjoint vertex lists (line_a, line_b),
 * this function returns only those edges that connect the sets. Each returned edge is
 * oriented as (a_label, b_label): if an original edge is (u in B, v in A), it is flipped
 * to (v, u).
 *
 * @tparam IntT Integral type for vertex labels.
 * @param global_edges The full set of undirected edges as pairs of original vertex labels.
 * @param line_a Vertex labels for layer A (original labels).
 * @param line_b Vertex labels for layer B (original labels).
 * @return std::vector<std::pair<IntT,IntT>> Interlayer edges oriented A→B (original labels).
 *
 * @note
 * - Edges with both endpoints in A or both in B are ignored.
 * - If the input contains both (u,v) and (v,u), both may appear (both normalized to A→B).
 *   See below for a dedup hint.
 * - Time complexity: O(E) average, where E = global_edges.size().
 */
  template <typename IntT>
  std::vector<std::pair<IntT,IntT>>
  ExtractInterlayerEdgesOriented(const std::vector<std::pair<IntT,IntT>>& global_edges,
				 const std::vector<IntT>& line_a,
				 const std::vector<IntT>& line_b) {
    std::unordered_set<IntT> A(line_a.begin(), line_a.end());
    std::unordered_set<IntT> B(line_b.begin(), line_b.end());
    
    std::vector<std::pair<IntT,IntT>> inter_edges;
    inter_edges.reserve(global_edges.size());
    
    for (const auto& [u, v] : global_edges) {
      const bool uA = A.count(u), vA = A.count(v);
      const bool uB = B.count(u), vB = B.count(v);
      
      // Keep only cross edges, and orient to (A, B)
      if (uA && vB) {
	inter_edges.emplace_back(u, v);        // already A→B
      } else if (uB && vA) {
	inter_edges.emplace_back(v, u);        // flip to A→B
      }
    }
    return inter_edges;
  }

/**
 * @brief Extract interlayer edges reindexed and oriented from A to B.
 *
 * Always returns pairs (a_idx, b_idx) where a_idx ∈ [0, |A|), b_idx ∈ [0, |B|).
 *
 * @tparam IntT Integral type for vertex labels.
 * @param global_edges The full set of edges in the graph (undirected), as pairs of vertex labels.
 * @param line_a Vertex labels for layer A.
 * @param line_b Vertex labels for layer B.
 * @return std::vector<std::pair<IntT,IntT>> Interlayer edges oriented A→B with local indices.
 *
 * @complexity O(E).
 */
  template <typename IntT>
  std::vector<std::pair<IntT,IntT>>
  ExtractInterlayerEdgesReindexedOriented(const std::vector<std::pair<IntT,IntT>>& global_edges,
					  const std::vector<IntT>& line_a,
					  const std::vector<IntT>& line_b) {
    std::unordered_map<IntT, IntT> a_local, b_local;
    a_local.reserve(line_a.size());
    b_local.reserve(line_b.size());
    for (size_t i = 0; i < line_a.size(); ++i) a_local[line_a[i]] = static_cast<IntT>(i);
    for (size_t i = 0; i < line_b.size(); ++i) b_local[line_b[i]] = static_cast<IntT>(i);
    
    std::vector<std::pair<IntT,IntT>> inter_edges;
    inter_edges.reserve(global_edges.size());
    
    for (const auto& [u,v] : global_edges) {
      auto iuA = a_local.find(u), ivA = a_local.find(v);
      auto iuB = b_local.find(u), ivB = b_local.find(v);
      
      if (iuA != a_local.end() && ivB != b_local.end()) {
	inter_edges.emplace_back(iuA->second, ivB->second);
      } else if (iuB != b_local.end() && ivA != a_local.end()) {
	inter_edges.emplace_back(ivA->second, iuB->second);
      }
    }
    return inter_edges;
  }

  /**
   * @brief Extract a subset of edges connected to a given set of sites.
   *
   * This function scans the input edge list `edges` and extracts all edges
   * that are incident to at least one vertex contained in `sites`.
   *
   * More precisely, for each edge (i,j) in `edges`, the edge is included
   * in the output if either i in sites or j in sites.
   *
   * The order of edges in the returned vector is preserved from the original
   * `edges` array. That is, if the extracted edges correspond to indices
   * m_1, m_2, ..., m_K in `edges`, then m_1 < m_2 < ... < m_K holds.
   *
   * @tparam IntT Integer type used for vertex labels.
   *
   * @param[in] edges
   *    List of edges represented as pairs (i,j), where i and j are vertex labels.
   *
   * @param[in] sites
   *    Subset of vertex labels. Edges connected to any of these vertices
   *    will be extracted.
   *
   * @return std::vector<std::pair<IntT,IntT>>
   *    A vector containing all edges from `edges` that are incident to at least
   *    one vertex in `sites`, preserving the original order.
   *
   * @note
   *    The membership test for `sites` is internally accelerated using
   *    std::unordered_set, resulting in average O(1) lookup time.
   *
   * @complexity
   *    Time complexity is O(|edges|+|sites|) on average.
   *
   * @warning
   *    This function selects edges where *at least one* endpoint is in `sites`.
   *    If you need edges where *both* endpoints are in `sites`, the condition
   *    should be modified accordingly.
   * 
   */
  template <typename IntT>
  std::vector<std::pair<IntT,IntT>>
  extract_incident_edges(const std::vector<std::pair<IntT,IntT>> & edges,
			 const std::vector<IntT> & sites) {
    std::unordered_set<IntT> site_set(sites.begin(),sites.end());
    std::vector<std::pair<IntT,IntT>> subedges;
    subedges.reserve(edges.size());
    for(const auto & e : edges) {
      if (site_set.count(e.first) || site_set.count(e.second)) {
	subedges.push_back(e);
      }
    }
    return subedges;
  }

  /**
   * @brief Extract edges whose both endpoints are contained in a given set of sites.
   *
   * This function scans the input edge list `edges` and extracts all edges
   * whose two endpoints are both contained in `sites`.
   *
   * More precisely, for each edge (i, j) in `edges`, the edge is included
   * in the output if both i ∈ sites and j ∈ sites.
   *
   * The order of edges in the returned vector is preserved from the original
   * `edges` array. That is, if the extracted edges correspond to indices
   * m₁, m₂, ..., m_K in `edges`, then m₁ < m₂ < ... < m_K holds.
   *
   * @tparam IntT Integer type used for vertex labels.
   *
   * @param[in] edges
   *   List of edges represented as pairs (i, j), where i and j are vertex labels.
   *
   * @param[in] sites
   *   Subset of vertex labels. Edges whose both endpoints are contained in
   *   this subset will be extracted.
   *
   * @return std::vector<std::pair<IntT, IntT>>
   *   A vector containing all edges from `edges` whose both endpoints are in
   *   `sites`, preserving the original order.
   *
   * @note
   *   This corresponds to extracting the edge set of the subgraph induced by
   *   `sites`, assuming `edges` represents the edge set of the original graph.
   *
   * @note
   *   The membership test for `sites` is internally accelerated using
   *   std::unordered_set, resulting in average O(1) lookup time.
   *
   * @complexity
   *   Time complexity is O(|edges| + |sites|) on average.
   *
   * @warning
   *   This function selects edges where *both* endpoints are in `sites`.
   *   If you need edges where *at least one* endpoint is in `sites`, use
   *   `extract_incident_edges` instead.
   */
  template <typename IntT>
  std::vector<std::pair<IntT, IntT>>
  extract_induced_edges(const std::vector<std::pair<IntT, IntT>>& edges,
			const std::vector<IntT>& sites) {
    std::unordered_set<IntT> site_set(sites.begin(), sites.end());
    
    std::vector<std::pair<IntT, IntT>> subedges;
    subedges.reserve(edges.size());
    
    for (const auto& e : edges) {
      if (site_set.count(e.first) && site_set.count(e.second)) {
	subedges.push_back(e);
      }
    }
    return subedges;
  }
  
  
}

#endif

