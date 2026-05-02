mpirun -np 1 ./boundary_mps --inputname input_graph.txt \
       --state_bond_dim 1 \
       --bmps_bp_tolerance 1.0e-8 --bmps_max_bp_steps 20 \
       --bmps_bond_dim 1 --bmps_tg_err 1.0e-8 --bmps_sv_min 1.0e-8 \
       --init_line 0,2,4,6 \
       --outputname output_graph_
