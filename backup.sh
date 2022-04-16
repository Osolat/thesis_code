#make clean
#make basic_bench
#make abe_circuit
#make gpsw
#make kLin_KP
cd objects
ulimit -s unlimited &&  ./gpsw_gap_ok 2 >> gpsw.txt && 
./gpsw_gap_ok 4 >> gpsw.txt && ./gpsw_gap_ok 8 >> gpsw.txt && ./gpsw_gap_ok 16 >> gpsw.txt && ./gpsw_gap_ok 32 >> gpsw.txt && ./gpsw_gap_ok 64 >> gpsw.txt && ./gpsw_gap_ok 128 >> gpsw.txt && ./gpsw_gap_ok 256 >> gpsw.txt && ./gpsw_gap_ok 512 >> gpsw.txt
#make clean
