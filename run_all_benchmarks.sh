#make clean
#make basic_bench
#make abe_circuit
#make gpsw
#make kLin_KP
cd objects
ulimit -s unlimited && 
./klin_kp_std 2 && ./klin_kp_G 2 && ./klin_kp_A 2 && ./klin_kp_P 2  && ./klin_kp_GAP 2 && ./klin_kp_GAP_OE 2 && ./klin_kp_GAP_OK 2 && 
./klin_kp_GAP_OD 2 && ./abe_circuit_base n_attr 2 && ./abe_circuit_gopt n_attr 2 && ./abe_circuit_api n_attr 2 &&   
./abe_circuit_pre n_attr 2 && ./abe_circuit_GAP n_attr 2 && ./abe_circuit_oe n_attr 2 && ./abe_circuit_ok n_attr 2 &&
./gpsw_ok 2 && ./gpsw_g_ok 2 && ./gpsw_a_ok 2 && ./gpsw_a_ok 2 && ./gpsw_gap_ok 2 && ./gpsw_gap_oe 2 ./gpsw_gap_od 2 &&
./klin_kp_std 32 && ./klin_kp_G 32 && ./klin_kp_A 32 && ./klin_kp_P 32  && ./klin_kp_GAP 32 && ./klin_kp_GAP_OE 32 && ./klin_kp_GAP_OK 32 && 
./klin_kp_GAP_OD 32 && ./abe_circuit_base n_attr 32 && ./abe_circuit_gopt n_attr 32 && ./abe_circuit_api n_attr 32 &&   
./abe_circuit_pre n_attr 32 && ./abe_circuit_GAP n_attr 32 && ./abe_circuit_oe n_attr 32 && ./abe_circuit_ok n_attr 32 &&
./gpsw_ok 32 && ./gpsw_g_ok 32 && ./gpsw_a_ok 32 && ./gpsw_p_ok 32 && ./gpsw_gap_ok 32 && ./gpsw_gap_oe 32 ./gpsw_gap_od 32 &&
./klin_kp_std 512 && ./klin_kp_G 512 && ./klin_kp_A 512 && ./klin_kp_P 512  && ./klin_kp_GAP 512 && ./klin_kp_GAP_OE 512 && ./klin_kp_GAP_OK 512 && 
./klin_kp_GAP_OD 512 && ./abe_circuit_base n_attr 512 && ./abe_circuit_gopt n_attr 512 && ./abe_circuit_api n_attr 512 &&   
./abe_circuit_pre n_attr 512 && ./abe_circuit_GAP n_attr 512 && ./abe_circuit_oe n_attr 512 && ./abe_circuit_ok n_attr 512 &&
./gpsw_ok 512 && ./gpsw_a_ok 512 && ./gpsw_a_ok 512 && ./gpsw_a_ok 512 && ./gpsw_gap_ok 512 && ./gpsw_gap_oe 512 ./gpsw_gap_od 512 &&
./klin_kp_std 1024 && ./klin_kp_G 1024 && ./klin_kp_A 1024 && ./klin_kp_P 1024  && ./klin_kp_GAP 1024 && ./klin_kp_GAP_OE 1024 && ./klin_kp_GAP_OK 1024 && 
./klin_kp_GAP_OD 1024 && ./abe_circuit_base n_attr 1024 && ./abe_circuit_gopt n_attr 1024 && ./abe_circuit_api n_attr 1024 &&   
./abe_circuit_pre n_attr 1024 && ./abe_circuit_GAP n_attr 1024 && ./abe_circuit_oe n_attr 1024 && ./abe_circuit_ok n_attr 1024 &&
./gpsw_ok 1024 && ./gpsw_g_ok 1024 && ./gpsw_g_ok 1024 && ./gpsw_p_ok 1024 && ./gpsw_gap_ok 1024 && ./gpsw_gap_oe 1024 ./gpsw_gap_od 1024
./bench_add_operations 100 && ./bench_iter_map_vs_map_sim 100 && ./bench_mul_opt 100 && ./bench_mul_pre 100 && ./bench_negations 100 && ./bench_operations 100 &&
./bench_sim_vs_pre 100 && ./bench_split_sim 50 &&
#make clean