/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/ptns/function/boundary.h
@brief Functions for boundary MPS contraction
*/
#ifndef TNBP_PTNS_FUNCTION_BOUNDARY_H
#define TBNP_PTNS_FUNCTION_BOUNDARY_H

#include <vector>
#include <unordered_set>
#include <utility> // for std::pair

#include "tnbp/framework/graph.h"

namespace tnbp {

  void BoundaryContractionLabelSetup(
       const std::vector<std::pair<int,int>> & global_edges,
       const std::vector<int> & initial_layer,
       std::vector<std::vector<int>> & layers,
       std::vector<std::vector<std::pair<int,int>>> & layer_edges) {

    layers = BuildContractLayers(global_edges,initial_layer);

    layer_edges.resize(layers.size());
    size_t m=0;
    for(const auto & line : layers) {
      layer_edges[m] = ExtractLineEdges(global_edges,line);
    }
  }

  template <typename TenT>
  void BoundaryContractionInitTensor(
       context_handle_t<TenT> & ctx,
       const std::vector<std::pair<int,int>> & global_edges,
       const std::vector<std::pair<int,int>> & boundary_edges,
       const std::vector<TenT> & V,
       const std::vector<int> & SiteIdx,
       const std::map<int,int> & Site_To_MpiRank,
       MPI_Comm comm,
       std::vector<TenT> & W) {
    
  }
  
}

#endif
