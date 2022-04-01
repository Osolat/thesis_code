
# Compilation makefile for all curves, all schemes
.PHONY: all clean

CXXFLAGS = -O3 -funroll-loops -fomit-frame-pointer -finline-small-functions -march=native -mtune=native 
CXX = g++

RELIC_LIB_BLS_12_381=/usr/local/lib/librelic_s.a
RELIC_INCLUDE_BLS_12_381=/usr/local/include/
#RELIC_LIB_BLS_12_381=~/Master-Thesis/relic/relic-target/lib/librelic_s.a
#RELIC_INCLUDE_BLS_12_381=~/Master-Thesis/

ARITH_OBJ = zp_arith.o g1_arith.o g2_arith.o gt_arith.o structures.o pairing_arith.o
LEGACY_OBJ = l_zobject.o l_zfunctioninput.o l_zattributelist.o l_zpolicy.o l_zdriver.o zscanner.o zparser.tab.o l_zgroup.o l_zelement.o l_zelement_bp.o l_zlsss.o
POLICY_OBJ = policy_tree.o
UTIL_OBJ = k_lin_util.o

main:
	$(info ************  Compiling ************)
	g++ $(CXXFLAGS) -c lib/k_lin/k_lin_util.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/policy/policy_tree.cpp -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/zp_arith.c -I$(RELIC_INCLUDE_BLS_12_381) 
	g++ $(CXXFLAGS) -c lib/g1_arith.c -I$(RELIC_INCLUDE_BLS_12_381) 
	g++ $(CXXFLAGS) -c lib/g2_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/gt_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/structures.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/pairing_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	$(CXX) -o objects/gpsw_g_ok $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) gpsw_g_ok.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	rm *.o

kLin_KP:
	$(info ************  Compiling K-LIN ************)
	g++ $(CXXFLAGS) -c lib/k_lin/k_lin_util.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/policy/policy_tree.cpp -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/zp_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/g1_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/g2_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/gt_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/structures.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/pairing_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	$(CXX) -o objects/klin_kp_1000_k2_std $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) kLin_KP.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/klin_kp_1000_k2_G $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) KLin_KP_G.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/klin_kp_1000_k2_A $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) KLin_KP_A.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/klin_kp_1000_k2_P $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) KLin_KP_P.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/klin_kp_1000_k2_GAP $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) KLin_KP_GAP.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/klin_kp_1000_k2_GAP_OE $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) kLin_KP_GAP_OE.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/klin_kp_1000_k2_GAP_OK $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) kLin_KP_GAP_OK.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/klin_kp_1000_k2_GAP_OD $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) kLin_KP_GAP_OD.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/bench_sim_vs_pre $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) bench_sim_vs_pre.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/bench_iter_map_vs_map_sim $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) bench_iter_map_vs_map_sim.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	rm *.o

abe_circuit:
	$(info ************  Compiling ************)
	g++ $(CXXFLAGS) -c lib/k_lin/k_lin_util.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/policy/policy_tree.cpp -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/zp_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/g1_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/g2_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/gt_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/structures.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c lib/pairing_arith.c -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c legacy/l_zobject.cpp
	g++ $(CXXFLAGS) -c legacy/l_zfunctioninput.cpp
	g++ $(CXXFLAGS) -c legacy/l_zattributelist.cpp
	g++ $(CXXFLAGS) -c legacy/l_zpolicy.cpp
	g++ $(CXXFLAGS) -c legacy/parser/l_zdriver.cpp
	g++ $(CXXFLAGS) -c legacy/parser/zscanner.cpp
	g++ $(CXXFLAGS) -c legacy/parser/zparser.tab.cc
	g++ $(CXXFLAGS) -c legacy/arith/l_zgroup.cpp
	cc $(CXXFLAGS) -I$(RELIC_INCLUDE_BLS_12_381) -fPIC  -O3 -DSSL_LIB_INIT  -Wno-implicit-function-declaration  -c legacy/arith/l_zelement.c -o l_zelement.o
	g++ $(CXXFLAGS) -c legacy/arith/l_zelement_bp.cpp -I$(RELIC_INCLUDE_BLS_12_381)
	g++ $(CXXFLAGS) -c legacy/lsss/l_zlsss.cpp -I$(RELIC_INCLUDE_BLS_12_381)
	$(CXX) -o objects/abe_circuit_oe $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(LEGACY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) abe_circuit_w_fanout_oe.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	rm *.o

policy_tree:
	$(info ************  Compiling Policy Tree Test************)
	g++ $(CXXFLAGS) -c lib/policy/policy_tree.cpp -I$(RELIC_INCLUDE_BLS_12_381)
	$(CXX) -o objects/main $(CXXFLAGS) $(POLICY_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) test_access_tree.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	rm *.o
	
benchmarks_ops:
	$(info ************  Benchmark for operations    ************)
	$(CXX) -o objects/main $(CXXFLAGS) -I$(RELIC_INCLUDE_BLS_12_381) bench_operations.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	rm *.o
failing_comp:
	$(info ************  Compiling Policy Tree Test************)
	$(CXX) -o objects/main $(CXXFLAGS) -I$(RELIC_INCLUDE_BLS_12_381) test_access_tree.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	rm *.o
clean:
	rm -f objects/* main *~

