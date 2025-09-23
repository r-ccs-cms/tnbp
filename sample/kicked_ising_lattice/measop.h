template <typename TenT>
void measops(tci::context_handle_t<TenT> & ctx,
	     const std::vector<std::pair<int,int>> & edges,
	     std::vector<int> & meas_sites,
	     std::vector<TenT> & Ox,
	     std::vector<TenT> & Oz,
	     std::vector<std::pair<int,int>> & meas_edges,
	     std::vector<TenT> & Oj) {

  using ElemT = typename tci::tensor_traits<TenT>::elem_t;
  using RealT = typename tci::tensor_traits<TenT>::real_t;
  using RankT = typename tci::tensor_traits<TenT>::rank_t;
  using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
  using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;

  auto sites = tnbp::GetSiteIndexFromBond(edges);
  
  ShapeT shape_z(2,2);
  std::vector<ElemT> data_z =
    { ElemT(1.0), ElemT(0.0),
      ElemT(0.0), ElemT(-1.0) };
  auto it_data_z = data_z.begin();
  TenT Zi = tci::assign_from_container<TenT>(
	       ctx,shape_z,it_data_z,
	       [](const CoorsT & coors) {
		 return coors[0]+coors[1]*2;
	       });

  ShapeT shape_x(2,2);
  std::vector<ElemT> data_x =
    { ElemT(0.0), ElemT(1.0),
      ElemT(1.0), ElemT(0.0) };
  auto it_data_x = data_x.begin();
  TenT Xi = tci::assign_from_container<TenT>(
	       ctx,shape_x,it_data_x,
	       [](const CoorsT & coors) {
		 return coors[0]+coors[1]*2;
	       });

  ShapeT shape_j(4,2);
  std::vector<ElemT> data_j =
    { ElemT(1.0), ElemT( 0.0), ElemT( 0.0), ElemT(0.0),
      ElemT(0.0), ElemT(-1.0), ElemT( 0.0), ElemT(0.0),
      ElemT(0.0), ElemT( 0.0), ElemT(-1.0), ElemT(0.0),
      ElemT(0.0), ElemT( 0.0), ElemT( 0.0), ElemT(1.0) };
  auto it_data_j = data_j.begin();
  TenT Ji = tci::assign_from_container<TenT>(
	       ctx,shape_j,it_data_j,
	       [](const CoorsT & coors) {
		 return coors[0]+coors[1]*2+coors[2]*4+coors[3]*8;
	       });

  Ox.resize(sites.size());
  Oz.resize(sites.size());
  Oj.resize(edges.size());

  for(auto & Oi : Ox) {
    tci::copy(ctx,Xi,Oi);
  }

  for(auto & Oi : Oz) {
    tci::copy(ctx,Zi,Oi);
  }

  for(auto & Oi : Oj) {
    tci::copy(ctx,Ji,Oi);
  }

  meas_sites = sites;
  meas_edges = edges;

}
			  
