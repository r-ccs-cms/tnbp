mpirun -np 1 ./estimator --backend ibm_kobe --circuit kicked_ising.qasm --max_bp_iterations 50 --sparse_pauli H_ising_kobe.txt --sv_min 1.0e-4 --truncation_error 1.0e-6
