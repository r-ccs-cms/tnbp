// coor -> address (row-major)
template <typename ShapeT, typename CoorT>
inline std::ptrdiff_t address_from_coor(const ShapeT& shape,
					const CoorT& coor) {
  if (shape.empty()) return 0;
  size_t addr = coor[shape.size()-1];
  for (size_t k = shape.size()-1; k > 0; --k) {
    addr = addr * shape[k-1] + coor[k-1];
  }
  return static_cast<std::ptrdiff_t>(addr);
}


template <typename TenT>
std::vector<TenT> CircuitTPO(
     tcapi::context_handle_t<TenT> & ctx,
     const std::vector<std::pair<int,int>> & edges,
     tcapi::real_t<TenT> Jz,
     tcapi::real_t<TenT> hz,
     tcapi::real_t<TenT> hx,
     tcapi::real_t<TenT> dt,
     MPI_Comm comm) {

  using ElemT = typename tcapi::tensor_traits<TenT>::elem_t;
  using RealT = typename tcapi::tensor_traits<TenT>::real_t;
  using OrderT = typename tcapi::tensor_traits<TenT>::order_t;
  using BondDimT = typename tcapi::tensor_traits<TenT>::bond_dim_t;
  using BondIdxT = typename tcapi::tensor_traits<TenT>::bond_idx_t;
  using BondLabelT = typename tcapi::tensor_traits<TenT>::bond_label_t;
  using ShapeT = typename tcapi::tensor_traits<TenT>::shape_t;
  using CoorsT = typename tcapi::tensor_traits<TenT>::elem_coors_t;
  using RealTenT = typename tcapi::tensor_traits<TenT>::real_ten_t;

  using ContextHandleR = typename tcapi::tensor_traits<RealTenT>::context_handle_t;
  ContextHandleR ctx_r;
  tcapi::create_context(ctx_r);

  int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
  int mpi_size; MPI_Comm_size(comm,&mpi_size);

  auto sites = tnbp::GetSiteIndexFromBond(edges);

  ShapeT shape_j(4,2);
  std::vector<ElemT> data_j =
    { ElemT(cos(dt*Jz),-sin(dt*Jz)), ElemT(0.0), ElemT(0.0), ElemT(0.0),
      ElemT(0.0), ElemT(cos(dt*Jz), sin(dt*Jz)), ElemT(0.0), ElemT(0.0),
      ElemT(0.0), ElemT(0.0), ElemT(cos(dt*Jz), sin(dt*Jz)), ElemT(0.0),
      ElemT(0.0), ElemT(0.0), ElemT(0.0), ElemT(cos(dt*Jz),-sin(dt*Jz)) };
  auto it_data_j = data_j.begin();
  TenT Uj = tcapi::assign_from_range<TenT>(
	      ctx,shape_j,it_data_j,
	      [](const CoorsT & coors) {
		return coors[0]+coors[1]*2+coors[2]*4+coors[3]*8;
	      });
  // Decomposed it to site tensors
  tcapi::List<BondIdxT> new_order_j = { 0, 2, 1, 3 };
  tcapi::transpose(ctx,Uj,new_order_j);
  TenT Ua;
  TenT Ub;
  RealTenT S;
  OrderT lb = 2;
  tcapi::svd(ctx,Uj,lb,Ua,S,Ub);
  tcapi::for_each(ctx_r,S,[](auto & elem) { elem = std::sqrt(elem); } );
  TenT D;
  if constexpr (std::is_same_v<TenT,RealTenT>) {
    D = tcapi::copy(ctx_r,S);
  } else {
    D = tcapi::to_cplx(ctx_r,S);
  }
  tcapi::List<BondLabelT> label_ua = {0,1,-1};
  tcapi::List<BondLabelT> label_da = {-1,2};
  tcapi::List<BondLabelT> label_va = {2,0,1};
  TenT Va;
  tcapi::contract(ctx,Ua,label_ua,D,label_da,Va,label_va);
  tcapi::List<BondLabelT> label_db = {0,-1};
  tcapi::List<BondLabelT> label_ub = {-1,1,2};
  tcapi::List<BondLabelT> label_vb = {0,1,2};
  TenT Vb;
  tcapi::contract(ctx,D,label_db,Ub,label_ub,Vb,label_vb);

  ShapeT shape_z(2,2);
  std::vector<ElemT> data_z =
    { ElemT(cos(dt*hz),-sin(dt*hz)), ElemT(0.0),
      ElemT(0.0), ElemT(cos(dt*hz), sin(dt*hz)) };
  auto it_data_z = data_z.begin();
  TenT Uz = tcapi::assign_from_range<TenT>(
	      ctx,shape_z,it_data_z,
	      [](const CoorsT & coors) {
		return coors[0]+coors[1]*2;
	      });

  ShapeT shape_x(2,2);
  std::vector<ElemT> data_x =
    { ElemT(cos(dt*hx)), ElemT(0.0,-sin(dt*hx)),
      ElemT(0.0,-sin(dt*hx)), ElemT(cos(dt*hx)) };
  auto it_data_x = data_x.begin();
  TenT Ux = tcapi::assign_from_range<TenT>(
	      ctx,shape_x,it_data_x,
	      [](const CoorsT & coors) {
		return coors[0]+coors[1]*2;
	      });

  // construct initial state
  
  std::vector<ElemT> data_i =
    { ElemT(1.0), ElemT(0.0),
      ElemT(0.0), ElemT(1.0) };

  std::vector<TenT> V(sites.size());

  for(size_t addr=0; addr < sites.size(); addr++) {
    auto bonds = tnbp::GetSurroundingBondIndex(sites[addr],edges);
    ShapeT shape_v(bonds.size()+2,1);
    shape_v[bonds.size()+0] = 2;
    shape_v[bonds.size()+1] = 2;
    auto it_data_i = data_i.begin();
    V[sites[addr]] = tcapi::assign_from_range<TenT>(
	     ctx,shape_v,it_data_i,
	     [&shape_v](const CoorsT & coors) {
	       return address_from_coor(shape_v,coors);
	     });
  }

  // j terms
  for(size_t m=0; m < edges.size(); m++) {
    TenT W;
    // a-site
    auto site_a = edges[m].first;
    auto order_ta = tcapi::order(ctx,V[site_a]);
    auto order_va = tcapi::order(ctx,Va);
    auto order_wa = order_ta + 1;
    tcapi::List<BondLabelT> label_ta(order_ta);
    tcapi::List<BondLabelT> label_va(order_va);
    tcapi::List<BondLabelT> label_wa(order_wa);
    ShapeT shape_ta = tcapi::shape(ctx,V[site_a]);
    ShapeT shape_va = tcapi::shape(ctx,Va);
    ShapeT shape_wa(order_ta);
    auto bonds_a = tnbp::GetSurroundingBondIndex(site_a,edges);
    BondLabelT contract_label = 0;
    for(size_t k=0; k < bonds_a.size(); k++) {
      if( bonds_a[k] == m ) {
	label_ta[k] = contract_label++;
	label_va[0] = contract_label++;
	shape_wa[k] = shape_ta[k]*shape_va[0];
      } else {
	label_ta[k] = contract_label++;
	shape_wa[k] = shape_ta[k];
      }
    }
    label_ta[order_ta-2] = contract_label++;
    label_ta[order_ta-1] = -1;
    label_va[1] = -1;
    label_va[2] = contract_label++;
    shape_wa[order_ta-2] = shape_ta[order_ta-2];
    shape_wa[order_ta-1] = shape_ta[order_ta-1];
    std::iota(label_wa.begin(),label_wa.end(),0);
    tcapi::contract(ctx,V[site_a],label_ta,Va,label_va,W,label_wa);
    tcapi::reshape(ctx,W,shape_wa,V[site_a]);

    // b-site
    auto site_b = edges[m].second;
    auto order_tb = tcapi::order(ctx,V[site_b]);
    auto order_vb = tcapi::order(ctx,Vb);
    auto order_wb = order_tb + 1;
    tcapi::List<BondLabelT> label_tb(order_tb);
    tcapi::List<BondLabelT> label_vb(order_vb);
    tcapi::List<BondLabelT> label_wb(order_wb);
    ShapeT shape_tb = tcapi::shape(ctx,V[site_b]);
    ShapeT shape_vb = tcapi::shape(ctx,Vb);
    ShapeT shape_wb(order_tb);
    auto bonds_b = tnbp::GetSurroundingBondIndex(site_b,edges);
    contract_label = 0;
    for(size_t k=0; k < bonds_b.size(); k++) {
      if( bonds_b[k] == m ) {
	label_tb[k] = contract_label++;
	label_vb[0] = contract_label++;
	shape_wb[k] = shape_tb[k]*shape_vb[0];
      } else {
	label_tb[k] = contract_label++;
	shape_wb[k] = shape_tb[k];
      }
    }
    label_tb[order_tb-2] = contract_label++;
    label_tb[order_tb-1] = -1;
    label_vb[1] = -1;
    label_vb[2] = contract_label++;
    shape_wb[order_tb-2] = shape_tb[order_tb-2];
    shape_wb[order_tb-1] = shape_tb[order_tb-1];
    std::iota(label_wb.begin(),label_wb.end(),0);
    tcapi::contract(ctx,V[site_b],label_tb,Vb,label_vb,W,label_wb);
    tcapi::reshape(ctx,W,shape_wb,V[site_b]);
  }

  // z terms
  for(size_t i=0; i < sites.size(); i++) {
    auto order_v = tcapi::order(ctx,V[i]);
    auto order_z = tcapi::order(ctx,Uz);
    auto order_w = tcapi::order(ctx,V[i]);
    tcapi::List<BondLabelT> label_v(order_v);
    tcapi::List<BondLabelT> label_z(order_z);
    tcapi::List<BondLabelT> label_w(order_v);
    std::iota(label_v.begin(),label_v.end(),0);
    label_v[order_v-1] = -1;
    label_z[0] = -1;
    label_z[1] = order_v-1;
    std::iota(label_w.begin(),label_w.end(),0);
    tcapi::contract(ctx,V[i],label_v,Uz,label_z,V[i],label_w);
  }

  // x terms
  for(size_t i=0; i < sites.size(); i++) {
    auto order_v = tcapi::order(ctx,V[i]);
    auto order_x = tcapi::order(ctx,Ux);
    auto order_w = tcapi::order(ctx,V[i]);
    tcapi::List<BondLabelT> label_v(order_v);
    tcapi::List<BondLabelT> label_x(order_x);
    tcapi::List<BondLabelT> label_w(order_w);
    std::iota(label_v.begin(),label_v.end(),0);
    label_v[order_v-1] = -1;
    label_x[0] = -1;
    label_x[1] = order_v-1;
    std::iota(label_w.begin(),label_w.end(),0);
    tcapi::contract(ctx,V[i],label_v,Ux,label_x,V[i],label_w);
  }

  return V;
  
}
					 
