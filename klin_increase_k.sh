#make clean
#make basic_bench
#make abe_circuit
#make gpsw
#make kLin_KP
cd objects
ulimit -s unlimited &&
#./klin_kp_std_k1 64 >> klin_k_inc.txt && ./klin_kp_G_k1 64 >> klin_k_inc.txt && ./klin_kp_A_k1 64 >> klin_k_inc.txt && 
#./klin_kp_GAP_OE_k1 64 >> klin_k_inc.txt && ./klin_kp_GAP_OK_k1 64 >> klin_k_inc.txt && 
#./klin_kp_GAP_OD_k1 64 >> klin_k_inc.txt &&
./klin_kp_std_2 64 >> klin_k_inc.txt && ./klin_kp_G_2 64 >> klin_k_inc.txt && ./klin_kp_A_2 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OE_2 64 >> klin_k_inc.txt && ./klin_kp_GAP_OK_2 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OD_2 64 >> klin_k_inc.txt &&
./klin_kp_std_3 64 >> klin_k_inc.txt && ./klin_kp_G_3 64 >> klin_k_inc.txt && ./klin_kp_A_3 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OE_3 64 >> klin_k_inc.txt && ./klin_kp_GAP_OK_3 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OD_3 64 >> klin_k_inc.txt &&
./klin_kp_std_4 64 >> klin_k_inc.txt && ./klin_kp_G_4 64 >> klin_k_inc.txt && ./klin_kp_A_4 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OE_4 64 >> klin_k_inc.txt && ./klin_kp_GAP_OK_4 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OD_4 64 >> klin_k_inc.txt &&
./klin_kp_std_5 64 >> klin_k_inc.txt && ./klin_kp_G_5 64 >> klin_k_inc.txt && ./klin_kp_A_5 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OE_5 64 >> klin_k_inc.txt && ./klin_kp_GAP_OK_5 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OD_5 64 >> klin_k_inc.txt &&
./klin_kp_std_6 64 >> klin_k_inc.txt && ./klin_kp_G_6 64 >> klin_k_inc.txt && ./klin_kp_A_6 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OE_6 64 >> klin_k_inc.txt && ./klin_kp_GAP_OK_6 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OD_6 64 >> klin_k_inc.txt &&
./klin_kp_std_7 64 >> klin_k_inc.txt && ./klin_kp_G_7 64 >> klin_k_inc.txt && ./klin_kp_A_7 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OE_7 64 >> klin_k_inc.txt && ./klin_kp_GAP_OK_7 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OD_7 64 >> klin_k_inc.txt &&
./klin_kp_std_8 64 >> klin_k_inc.txt && ./klin_kp_G_8 64 >> klin_k_inc.txt && ./klin_kp_A_8 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OE_8 64 >> klin_k_inc.txt && ./klin_kp_GAP_OK_8 64 >> klin_k_inc.txt && 
./klin_kp_GAP_OD_8 64 >> klin_k_inc.txt &&
#./klin_kp_std_k9 64 >> klin_k_inc.txt && ./klin_kp_G_k9 64 >> klin_k_inc.txt && ./klin_kp_A_k9 64 >> klin_k_inc.txt && 
#./klin_kp_GAP_OE_k9 64 >> klin_k_inc.txt && ./klin_kp_GAP_OK_k9 64 >> klin_k_inc.txt && 
#./klin_kp_GAP_OD_k9 64 >> klin_k_inc.txt &&
#./klin_kp_std_k10 64 >> klin_k_inc.txt && ./klin_kp_G_k10 64 >> klin_k_inc.txt && ./klin_kp_A_k10 64 >> klin_k_inc.txt && 
#./klin_kp_GAP_OE_k10 64 >> klin_k_inc.txt && ./klin_kp_GAP_OK_k10 64 >> klin_k_inc.txt && 
#./klin_kp_GAP_OD_k10 64 >> klin_k_inc.txt
./klin_kp_GAP_OE_2 1024 >> klin_k_inc.txt && ./klin_kp_GAP_OK_2 1024 >> klin_k_inc.txt && 
./klin_kp_GAP_OD_2 1024 >> klin_k_inc.txt 


#make clean
