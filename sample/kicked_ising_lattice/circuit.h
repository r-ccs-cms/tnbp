template <typename TenT>
std::vector<TenT> CircuitTPO(
     tci::context_handle_t<TenT> & ctx,
     int num_qubits,
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
  tci::create_context<RealTenT>(ctx_r);

  int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
  int mpi_size; MPI_Comm_size(comm,&mpi_size);

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

  // construct initial TPO
  ShapeT shape_e(3,2);
  shape_e[0] = 1;
  ShapeT shape_i(4,2);
  shape_i[0] = 1;
  shape_i[1] = 1;
  std::vector<ElemT> data_i =
    { ElemT(1.0), ElemT(0.0),
      ElemT(0.0), ElemT(1.0) };
  auto it_data_i = data_i.begin();
  TenT Ve = tci::assign_from_container<TenT>(
	       ctx,shape_e,it_data_i,
	       [](const CoorsT & coors) {
		 return coors[0]+coors[1]*1+coors[2]*2;
	       } );
  TenT Vi = tci::assign_from_container<TenT>(
	       ctx,shape_i,it_data_i,
	       [](const CoorsT & coors) {
		 return coors[0]+coors[1]*1+coors[2]*1+coors[3]*2;
	       } );

  std::vector<TenT> V(num_qubits);
  for(int i=0; i < num_qubits; i++) {
    if( i == 0 || i == num_qubits-1 ) {
      tci::copy(ctx,Ve,V[i]);
    } else {
      tci::copy(ctx,Vi,V[i]);
    }
  }

  // x terms
  
  for(int i=0; i < num_qubits; i++) {
    if( i == 0 || i == num_qubits-1 ) {
      tci::List<BondLabelT> label_v = {0,1,-1};
      tci::List<BondLabelT> label_x = {-1,2};
      tci::List<BondLabelT> label_w = {0,1,2};
      tci::contract(ctx,V[i],label_v,Ux,label_x,V[i],label_w);
    } else {
      tci::List<BondLabelT> label_v = {0,1,2,-1};
      tci::List<BondLabelT> label_x = {-1,3};
      tci::List<BondLabelT> label_w = {0,1,2,3};
      tci::contract(ctx,V[i],label_v,Ux,label_x,V[i],label_w);
    }
  }
  
  // j terms
  for(int i=0; i < num_qubits-1; i++) {
    TenT W;
    if ( i == 0 ) {
      tci::List<BondLabelT> label_va = {0,2,-1};
      tci::List<BondLabelT> label_a = {1,-1,3};
      tci::List<BondLabelT> label_wa = {0,1,2,3};
      ShapeT shape_va = tci::shape(ctx,V[i]);
      ShapeT shape_a = tci::shape(ctx,Va);
      ShapeT new_shape_va = {shape_va[0]*shape_a[0],shape_va[1],shape_va[2]};
      tci::contract(ctx,V[i],label_va,Va,label_a,W,label_wa);
      tci::reshape(ctx,W,new_shape_va,V[i]);
      tci::List<BondLabelT> label_vb = {0,2,3,-1};
      tci::List<BondLabelT> label_b = {1,-1,4};
      tci::List<BondLabelT> label_wb = {0,1,2,3,4};
      ShapeT shape_vb = tci::shape(ctx,V[i+1]);
      ShapeT shape_b = tci::shape(ctx,Vb);
      ShapeT new_shape_vb = {shape_vb[0]*shape_b[0],shape_vb[1],shape_vb[2],shape_vb[3]};
      tci::contract(ctx,V[i+1],label_vb,Vb,label_b,W,label_wb);
      tci::reshape(ctx,W,new_shape_vb,V[i+1]);
    } else if ( i == num_qubits-2 ) {
      tci::List<BondLabelT> label_va = {0,1,3,-1};
      tci::List<BondLabelT> label_a = {2,-1,4};
      tci::List<BondLabelT> label_wa = {0,1,2,3,4};
      ShapeT shape_va = tci::shape(ctx,V[i]);
      ShapeT shape_a = tci::shape(ctx,Va);
      ShapeT new_shape_va = {shape_va[0],shape_va[1]*shape_a[0],shape_va[2],shape_va[3]};
      tci::contract(ctx,V[i],label_va,Va,label_a,W,label_wa);
      tci::reshape(ctx,W,new_shape_va,V[i]);
      tci::List<BondLabelT> label_vb = {0,2,-1};
      tci::List<BondLabelT> label_b = {1,-1,3};
      tci::List<BondLabelT> label_wb = {0,1,2,3};
      ShapeT shape_vb = tci::shape(ctx,V[i+1]);
      ShapeT shape_b = tci::shape(ctx,Vb);
      ShapeT new_shape_vb = {shape_vb[0]*shape_b[0],shape_vb[1],shape_vb[2]};
      tci::contract(ctx,V[i+1],label_vb,Vb,label_b,W,label_wb);
      tci::reshape(ctx,W,new_shape_vb,V[i+1]);
    } else {
      tci::List<BondLabelT> label_va = {0,1,3,-1};
      tci::List<BondLabelT> label_a = {2,-1,4};
      tci::List<BondLabelT> label_wa = {0,1,2,3,4};
      ShapeT shape_va = tci::shape(ctx,V[i]);
      ShapeT shape_a = tci::shape(ctx,Va);
      ShapeT new_shape_va = {shape_va[0],shape_va[1]*shape_a[0],shape_va[2],shape_va[3]};
      tci::contract(ctx,V[i],label_va,Va,label_a,W,label_wa);
      tci::reshape(ctx,W,new_shape_va,V[i]);
      tci::List<BondLabelT> label_vb = {0,2,3,-1};
      tci::List<BondLabelT> label_b = {1,-1,4};
      tci::List<BondLabelT> label_wb = {0,1,2,3,4};
      ShapeT shape_vb = tci::shape(ctx,V[i+1]);
      ShapeT shape_b = tci::shape(ctx,Vb);
      ShapeT new_shape_vb = {shape_vb[0]*shape_b[0],shape_vb[1],shape_vb[2],shape_vb[3]};
      tci::contract(ctx,V[i+1],label_vb,Vb,label_b,W,label_wb);
      tci::reshape(ctx,W,new_shape_vb,V[i+1]);
    }
  }

  // z terms
  for(int i=0; i < num_qubits; i++) {
    if( i == 0 || i == num_qubits-1 ) {
      tci::List<BondLabelT> label_v = {0,1,-1};
      tci::List<BondLabelT> label_z = {-1,2};
      tci::List<BondLabelT> label_w = {0,1,2};
      tci::contract(ctx,V[i],label_v,Uz,label_z,V[i],label_w);
    } else {
      tci::List<BondLabelT> label_v = {0,1,2,-1};
      tci::List<BondLabelT> label_z = {-1,3};
      tci::List<BondLabelT> label_w = {0,1,2,3};
      tci::contract(ctx,V[i],label_v,Uz,label_z,V[i],label_w);
    }
  }

  return V;
  
}
					 
