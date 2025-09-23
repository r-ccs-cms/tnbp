mpirun -np 4 ./estimator --backend ibm_kobe --circuit kicked_ising.qasm --num_gates 840 --max_bp_iterations 50 --sparse_pauli H_ising_kobe.txt
