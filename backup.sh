#make clean
#make basic_bench
#make abe_circuit
#make gpsw
#make kLin_KP
cd objects
ulimit -s unlimited &&  
./klin_kp_A 2 >> klin.txt && ./klin_kp_A 4 >> klin.txt && ./klin_kp_A 8 >> klin.txt && ./klin_kp_A 16 >> klin.txt && ./klin_kp_A 32 >> klin.txt && ./klin_kp_A 64 >> klin.txt && ./klin_kp_A 128 >> klin.txt && ./klin_kp_A 256 >> klin.txt && ./klin_kp_A 512 >> klin.txt && ./klin_kp_P 2 >> klin.txt && ./klin_kp_P 4 >> klin.txt && 
./klin_kp_GAP 2 >> klin.txt && ./klin_kp_GAP 4 >> klin.txt && ./klin_kp_GAP 8 >> klin.txt && ./klin_kp_GAP 16 >> klin.txt && 
./klin_kp_GAP 32 >> klin.txt && ./klin_kp_GAP 64 >> klin.txt && ./klin_kp_GAP 128 >> klin.txt && ./klin_kp_GAP 256 >> klin.txt && ./klin_kp_GAP 512 >> klin.txt && ./klin_kp_GAP_OE 2 >> klin.txt && ./klin_kp_GAP_OE 4 >> klin.txt && ./klin_kp_GAP_OE 8 >> klin.txt && ./klin_kp_GAP_OE 16 >> klin.txt && ./klin_kp_GAP_OE 32 >> klin.txt && 
./klin_kp_GAP_OE 64 >> klin.txt && ./klin_kp_GAP_OE 128 >> klin.txt && ./klin_kp_GAP_OE 256 >> klin.txt && ./klin_kp_GAP_OE 512 >> klin.txt && ./klin_kp_GAP_OK 2 >> klin.txt && ./klin_kp_GAP_OK 4 >> klin.txt && ./klin_kp_GAP_OK 8 >> klin.txt && ./klin_kp_GAP_OK 16 >> klin.txt && ./klin_kp_GAP_OK 32 >> klin.txt && ./klin_kp_GAP_OK 64 >> klin.txt && 
./klin_kp_GAP_OK 128 >> klin.txt && ./klin_kp_GAP_OK 256 >> klin.txt && ./klin_kp_GAP_OK 512 >> klin.txt && ./klin_kp_GAP_OD 2 >> klin.txt && ./klin_kp_GAP_OD 4 >> klin.txt && ./klin_kp_GAP_OD 8 >> klin.txt && ./klin_kp_GAP_OD 16 >> klin.txt && ./klin_kp_GAP_OD 32 >> klin.txt && ./klin_kp_GAP_OD 64 >> klin.txt && ./klin_kp_GAP_OD 128 >> klin.txt && 
./klin_kp_GAP_OD 256 >> klin.txt && ./klin_kp_GAP_OD 512 >> klin.txt && 
./bench_iter_map_vs_map_sim  256 >> basic_bench.txt && 
./bench_mul_opt  256 >> basic_bench.txt && ./bench_mul_pre  256 >> basic_bench.txt && 
./bench_negations  256 >> basic_bench.txt &&
./klin_kp_A 1024 >> klin.txt && ./klin_kp_GAP 1024 >> klin.txt && ./klin_kp_GAP_OE 1024 >> klin.txt && ./klin_kp_GAP_OK 1024 >> klin.txt && ./klin_kp_GAP_OD 1024 >> klin.txt
#make clean
