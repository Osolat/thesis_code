#make clean
#make basic_bench
#make abe_circuit
#make gpsw
#make kLin_KP
cd objects
ulimit -s unlimited &&
./klin_kp_GAP_OD_2 1024 >> klin_k_inc.txt && ./klin_kp_GAP_OE_2 1024 >> klin_k_inc.txt && ./klin_kp_GAP_OK_2 1024 >> klin_k_inc.txt

#make clean
