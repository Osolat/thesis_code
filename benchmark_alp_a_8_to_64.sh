cd objects
ulimit -s unlimited &&
./bench_simultaneous_stuff >> sim_bench.txt &&
#./alp_naive n_attr 2 >> alp_result.txt &&
#./alp_g_oe n_attr 2 >> alp_result.txt &&
#./alp_a_oe n_attr 2 >> alp_result.txt &&
./alp_p_oe n_attr 2 >> alp_result.txt &&
./alp_GAP_oe n_attr 2 >> alp_result.txt &&
#./alp_GAP_ok n_attr 2 >> alp_result.txt &&
#./alp_naive n_attr 4 >> alp_result.txt &&
#./alp_g_oe n_attr 4 >> alp_result.txt &&
#./alp_a_oe n_attr 4 >> alp_result.txt &&
./alp_p_oe n_attr 4 >> alp_result.txt &&
./alp_GAP_oe n_attr 4 >> alp_result.txt &&
#./alp_GAP_ok n_attr 4 >> alp_result.txt &&
#./alp_naive n_attr 8 >> alp_result.txt &&
#./alp_g_oe n_attr 8 >> alp_result.txt &&
#./alp_a_oe n_attr 8 >> alp_result.txt &&
./alp_p_oe n_attr 8 >> alp_result.txt &&
./alp_GAP_oe n_attr 8 >> alp_result.txt &&
#./alp_GAP_ok n_attr 8 >> alp_result.txt &&
#./alp_naive n_attr 16 >> alp_result.txt &&
#./alp_g_oe n_attr 16 >> alp_result.txt &&
#./alp_a_oe n_attr 16 >> alp_result.txt &&
./alp_p_oe n_attr 16 >> alp_result.txt &&
./alp_GAP_oe n_attr 16 >> alp_result.txt &&
#./alp_GAP_ok n_attr 16 >> alp_result.txt &&
#./alp_naive n_attr 32 >> alp_result.txt &&
#./alp_g_oe n_attr 32 >> alp_result.txt &&
#./alp_a_oe n_attr 32 >> alp_result.txt &&
./alp_p_oe n_attr 32 >> alp_result.txt &&
./alp_GAP_oe n_attr 32 >> alp_result.txt &&
#./alp_GAP_ok n_attr 32 >> alp_result.txt &&
#./alp_naive n_attr 64 >> alp_result.txt &&
#./alp_g_oe n_attr 64 >> alp_result.txt &&
#./alp_a_oe n_attr 64 >> alp_result.txt &&
./alp_p_oe n_attr 64 >> alp_result.txt &&
./alp_GAP_oe n_attr 64 >> alp_result.txt &&
#./alp_GAP_ok n_attr 64 >> alp_result.txt 
./alp_p_oe_500 n_attr 128 >> alp_result.txt &&
./alp_GAP_oe_500 n_attr 128 >> alp_result.txt &&
./alp_p_oe_100 n_attr 256 >> alp_result.txt &&
./alp_GAP_oe_100 n_attr 256 >> alp_result.txt &&
./alp_p_oe_10 n_attr 512 >> alp_result.txt &&
./alp_GAP_oe_10 n_attr 512 >> alp_result.txt &&
./alp_p_oe_2 n_attr 1024 >> alp_result.txt &&
./alp_GAP_oe_2 n_attr 1024 >> alp_result.txt 