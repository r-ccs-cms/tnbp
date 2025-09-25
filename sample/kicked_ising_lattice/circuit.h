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
     tci::context_handle_t<TenT> & ctx,
     const std::vector<std::pair<int,int>> & edges,
     tci::real_t<TenT> Jz,
     tci::real_t<TenT> hz,
     tci::real_t<TenT> hx,
     tci::real_t<TenT> dt,
     MPI_Comm comm) {

  using ElemT = typename tci::tensor_traits<TenT>::elem_t;
  using RealT = typename tci::tensor_traits<TenT>::real_t;
  using RankT = typename tci::tensor_traits<TenT>::rank_t;
  using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;
  using BondIdxT = typename tci::tensor_traits<TenT>::bond_idx_t;
  using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
  using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
  using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
  using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;

  using ContextHandleR = typename tci::tensor_traits<RealTenT>::context_handle_t;
  ContextHandleR ctx_r;
  tci::create_context(ctx_r);

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
  TenT Uj = tci::assign_from_container<TenT>(
	      ctx,shape_j,it_data_j,
	      [](const CoorsT & coors) {
		return coors[0]+coors[1]*2+coors[2]*4+coors[3]*8;
	      });
  // Decomposed it to site tensors
  tci::List<BondIdxT> new_order_j = { 0, 2, 1, 3 };
  tci::transpose(ctx,Uj,new_order_j);
  TenT Ua;
  TenT Ub;
  RealTenT S;
  RankT lb = 2;
  tci::svd(ctx,Uj,lb,Ua,S,Ub);
  tci::for_each(ctx_r,S,[](RealT & elem) { elem = std::sqrt(elem); } );
  TenT D;
  tci::convert(ctx_r,S,ctx,D);
  tci::diag(ctx,D);
  tci::List<BondLabelT> label_ua = {0,1,-1};
  tci::List<BondLabelT> label_da = {-1,2};
  tci::List<BondLabelT> label_va = {2,0,1};
  TenT Va;
  tci::contract(ctx,Ua,label_ua,D,label_da,Va,label_va);
  tci::List<BondLabelT> label_db = {0,-1};
  tci::List<BondLabelT> label_ub = {-1,1,2};
  tci::List<BondLabelT> label_vb = {0,1,2};
  TenT Vb;
  tci::contract(ctx,D,label_db,Ub,label_ub,Vb,label_vb);

  ShapeT shape_z(2,2);
  std::vector<ElemT> data_z =
    { ElemT(cos(dt*hz),-sin(dt*hz)), ElemT(0.0),
      ElemT(0.0), ElemT(cos(dt*hz), sin(dt*hz)) };
  auto it_data_z = data_z.begin();
  TenT Uz = tci::assign_from_container<TenT>(
	      ctx,shape_z,it_data_z,
	      [](const CoorsT & coors) {
		return coors[0]+coors[1]*2;
	      });

  ShapeT shape_x(2,2);
  std::vector<ElemT> data_x =
    { ElemT(cos(dt*hx)), ElemT(0.0,-sin(dt*hx)),
      ElemT(0.0,-sin(dt*hx)), ElemT(cos(dt*hx)) };
  auto it_data_x = data_x.begin();
  TenT Ux = tci::assign_from_container<TenT>(
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
    V[sites[addr]] = tci::assign_from_container<TenT>(
	     ctx,shape_v,it_data_i,
	     [&shape_v](const CoorsT & coors) {
	       return address_from_coor(shape_v,coors);
	     });
  }

  std::cout << " End construction of initial state " << std::endl;
  
  // x terms
  
  for(size_t i=0; i < sites.size(); i++) {
    auto rank_v = tci::rank(ctx,V[i]);
    auto rank_x = tci::rank(ctx,Ux);
    auto rank_w = tci::rank(ctx,V[i]);
    tci::List<BondLabelT> label_v(rank_v);
    tci::List<BondLabelT> label_x(rank_x);
    tci::List<BondLabelT> label_w(rank_w);
    std::iota(label_v.begin(),label_v.end(),0);
    label_v[rank_v-1] = -1;
    label_x[0] = -1;
    label_x[1] = rank_v-1;
    std::iota(label_w.begin(),label_w.end(),0);
    tci::contract(ctx,V[i],label_v,Ux,label_x,V[i],label_w);
  }
  
  std::cout << " End construction of X-gates " << std::endl;
  
  // j terms
  for(size_t m=0; m < edges.size(); m++) {
    TenT W;
    // a-site
    auto site_a = edges[m].first;
    auto rank_ta = tci::rank(ctx,V[site_a]);
    auto rank_va = tci::rank(ctx,Va);
    auto rank_wa = rank_ta + 1;
    tci::List<BondLabelT> label_ta(rank_ta);
    tci::List<BondLabelT> label_va(rank_va);
    tci::List<BondLabelT> label_wa(rank_wa);
    ShapeT shape_ta = tci::shape(ctx,V[site_a]);
    ShapeT shape_va = tci::shape(ctx,Va);
    ShapeT shape_wa(rank_ta);
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
    label_ta[rank_ta-2] = contract_label++;
    label_ta[rank_ta-1] = -1;
    label_va[1] = -1;
    label_va[2] = contract_label++;
    shape_wa[rank_ta-2] = shape_ta[rank_ta-2];
    shape_wa[rank_ta-1] = shape_ta[rank_ta-1];
    std::iota(label_wa.begin(),label_wa.end(),0);
    tci::contract(ctx,V[site_a],label_ta,Va,label_va,W,label_wa);
    tci::reshape(ctx,W,shape_wa,V[site_a]);

    // b-site
    auto site_b = edges[m].second;
    auto rank_tb = tci::rank(ctx,V[site_b]);
    auto rank_vb = tci::rank(ctx,Vb);
    auto rank_wb = rank_tb + 1;
    tci::List<BondLabelT> label_tb(rank_tb);
    tci::List<BondLabelT> label_vb(rank_vb);
    tci::List<BondLabelT> label_wb(rank_wb);
    ShapeT shape_tb = tci::shape(ctx,V[site_b]);
    ShapeT shape_vb = tci::shape(ctx,Vb);
    ShapeT shape_wb(rank_tb);
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
    label_tb[rank_tb-2] = contract_label++;
    label_tb[rank_tb-1] = -1;
    label_vb[1] = -1;
    label_vb[2] = contract_label++;
    shape_wb[rank_tb-2] = shape_tb[rank_tb-2];
    shape_wb[rank_tb-1] = shape_tb[rank_tb-1];
    std::iota(label_wb.begin(),label_wb.end(),0);
    tci::contract(ctx,V[site_b],label_tb,Vb,label_vb,W,label_wb);
    tci::reshape(ctx,W,shape_wb,V[site_b]);
  }

  std::cout << " End construction of J-gates " << std::endl;
  
  // z terms
  for(size_t i=0; i < sites.size(); i++) {
    auto rank_v = tci::rank(ctx,V[i]);
    auto rank_z = tci::rank(ctx,Uz);
    auto rank_w = tci::rank(ctx,V[i]);
    tci::List<BondLabelT> label_v(rank_v);
    tci::List<BondLabelT> label_z(rank_z);
    tci::List<BondLabelT> label_w(rank_v);
    std::iota(label_v.begin(),label_v.end(),0);
    label_v[rank_v-1] = -1;
    label_z[0] = -1;
    label_z[1] = rank_v-1;
    std::iota(label_w.begin(),label_w.end(),0);
    tci::contract(ctx,V[i],label_v,Uz,label_z,V[i],label_w);
  }

  std::cout << " End construction of Z-gates " << std::endl;
  
  return V;
  
}
					 
