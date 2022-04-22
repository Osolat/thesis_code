cd objects
ulimit -s unlimited &&
./alp_GAP_oe n_attr 2 >> alp_result.txt &&
./alp_GAP_ok n_attr 2 >> alp_result.txt &&
./alp_GAP_oe n_attr 4 >> alp_result.txt &&
./alp_GAP_ok n_attr 4 >> alp_result.txt &&
./alp_GAP_oe n_attr 8 >> alp_result.txt &&
./alp_GAP_ok n_attr 8 >> alp_result.txt &&
./alp_GAP_oe n_attr 16 >> alp_result.txt &&
./alp_GAP_ok n_attr 16 >> alp_result.txt &&
./alp_GAP_oe n_attr 32 >> alp_result.txt &&
./alp_GAP_ok n_attr 32 >> alp_result.txt &&
./alp_GAP_oe n_attr 64 >> alp_result.txt &&
./alp_GAP_ok n_attr 64 >> alp_result.txt &&
./alp_GAP_oe n_attr 128 >> alp_result.txt &&
./alp_GAP_ok n_attr 128 >> alp_result.txt &&
./alp_GAP_oe n_attr 256 >> alp_result.txt &&
./alp_GAP_ok n_attr 256 >> alp_result.txt &&
./alp_GAP_oe n_attr 512 >> alp_result.txt &&
./alp_GAP_ok n_attr 512 >> alp_result.txt &&
./alp_GAP_oe n_attr 1024 >> alp_result.txt &&
./alp_GAP_ok n_attr 1024 >> alp_result.txt 