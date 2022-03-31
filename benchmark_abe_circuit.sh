make clean
make abe_circuit
ulimit -s unlimited && objects/abe_circuit_base n_attr 2 32 512 1024 && objects/abe_circuit_gopt n_attr 2 32 512 1024 && objects/abe_circuit_api n_attr 2 32 512 1024 &&   
objects/abe_circuit_pre n_attr 2 32 512 1024 && objects/abe_circuit_GAP n_attr 2 32 512 1024 && objects/abe_circuit_oe n_attr 2 32 512 1024 && objects/abe_circuit_ok
make clean