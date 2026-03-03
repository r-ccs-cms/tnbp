mpirun -np 1 ./qasm_to_tpo --circuit ising_circuit.qasm \
       --do_opt_tpo 0 --opt_tpo_eps 1.0e-8

