cd objects
ulimit -s unlimited &&
./alp_naive n_attr 8 >> alp_result.txt &&
./alp_g_oe n_attr 8 >> alp_result.txt &&
./alp_a_oe n_attr 8 >> alp_result.txt &&
./alp_p_oe n_attr 8 >> alp_result.txt &&
./alp_GAP_oe n_attr 8 >> alp_result.txt &&
./alp_GAP_ok n_attr 8 >> alp_result.txt &&
./alp_naive n_attr 16 >> alp_result.txt &&
./alp_g_oe n_attr 16 >> alp_result.txt &&
./alp_a_oe n_attr 16 >> alp_result.txt &&
./alp_p_oe n_attr 16 >> alp_result.txt &&
./alp_GAP_oe n_attr 16 >> alp_result.txt &&
./alp_GAP_ok n_attr 16 >> alp_result.txt &&
./alp_naive n_attr 32 >> alp_result.txt &&
./alp_g_oe n_attr 32 >> alp_result.txt &&
./alp_a_oe n_attr 32 >> alp_result.txt &&
./alp_p_oe n_attr 32 >> alp_result.txt &&
./alp_GAP_oe n_attr 32 >> alp_result.txt &&
./alp_GAP_ok n_attr 32 >> alp_result.txt &&
./alp_naive n_attr 64 >> alp_result.txt &&
./alp_g_oe n_attr 64 >> alp_result.txt &&
./alp_a_oe n_attr 64 >> alp_result.txt &&
./alp_p_oe n_attr 64 >> alp_result.txt &&
./alp_GAP_oe n_attr 64 >> alp_result.txt &&
./alp_GAP_ok n_attr 64 >> alp_result.txt 