mpirun -np 1 ./ptnos-bp --backend ibm_kobe --circuit kicked_ising.qasm \
       --measurement_barrier 3 \
       --sparse_pauli H_ising_kobe.txt \
       --max_bp_iterations 50 --sv_min 1.0e-4 --truncation_error 1.0e-6 \
       --step_start 0 --step_end 3 \
