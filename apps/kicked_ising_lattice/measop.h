template <typename TenT>
void measops(tcapi::context_handle_t<TenT> & ctx,
	     const std::vector<std::pair<int,int>> & edges,
	     std::vector<int> & meas_sites,
	     std::vector<TenT> & Ox,
	     std::vector<TenT> & Oz,
	     std::vector<std::pair<int,int>> & meas_edges,
	     std::vector<TenT> & Oj) {

  using ElemT = typename tcapi::tensor_traits<TenT>::elem_t;
  using RealT = typename tcapi::tensor_traits<TenT>::real_t;
  using OrderT = typename tcapi::tensor_traits<TenT>::order_t;
  using ShapeT = typename tcapi::tensor_traits<TenT>::shape_t;
  using CoorsT = typename tcapi::tensor_traits<TenT>::elem_coors_t;

  auto sites = tnbp::GetSiteIndexFromBond(edges);
  
  ShapeT shape_z(2,2);
  std::vector<ElemT> data_z =
    { ElemT(1.0), ElemT(0.0),
      ElemT(0.0), ElemT(-1.0) };
  auto it_data_z = data_z.begin();
  TenT Zi = tcapi::assign_from_range<TenT>(
	       ctx,shape_z,it_data_z,
	       [](const CoorsT & coors) {
		 return coors[0]+coors[1]*2;
	       });

  ShapeT shape_x(2,2);
  std::vector<ElemT> data_x =
    { ElemT(0.0), ElemT(1.0),
      ElemT(1.0), ElemT(0.0) };
  auto it_data_x = data_x.begin();
  TenT Xi = tcapi::assign_from_range<TenT>(
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
  TenT Ji = tcapi::assign_from_range<TenT>(
	       ctx,shape_j,it_data_j,
	       [](const CoorsT & coors) {
		 return coors[0]+coors[1]*2+coors[2]*4+coors[3]*8;
	       });

  Ox.resize(sites.size());
  Oz.resize(sites.size());
  Oj.resize(edges.size());

  for(auto & Oi : Ox) {
    Oi = tcapi::copy(ctx,Xi);
  }

  for(auto & Oi : Oz) {
    Oi = tcapi::copy(ctx,Zi);
  }

  for(auto & Oi : Oj) {
    Oi = tcapi::copy(ctx,Ji);
  }

  meas_sites = sites;
  meas_edges = edges;

}
			  
