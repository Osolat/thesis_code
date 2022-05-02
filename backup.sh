#make clean
#make basic_bench
#make abe_circuit
#make gpsw
#make kLin_KP
cd objects
ulimit -s unlimited && 
./bench_tree_coeffs_vs_rand 2 >> basic_bench.txt && ./bench_tree_coeffs_vs_rand 4 >> basic_bench.txt && ./bench_tree_coeffs_vs_rand 8 >> basic_bench.txt && ./bench_tree_coeffs_vs_rand 16 >> basic_bench.txt &&
./bench_tree_coeffs_vs_rand 32 >> basic_bench.txt && ./bench_tree_coeffs_vs_rand 64 >> basic_bench.txt && ./bench_tree_coeffs_vs_rand 128 >> basic_bench.txt && ./bench_tree_coeffs_vs_rand 256 >> basic_bench.txt &&
./bench_tree_coeffs_vs_rand 512 >> basic_bench.txt && ./bench_tree_coeffs_vs_rand 1024 >> basic_bench.txt &&
./gpsw_gap_ok 1024 >> gpsw.txt
#make clean
