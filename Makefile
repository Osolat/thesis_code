
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
	$(CXX) -o objects/gpsw_oe_10000 $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) run_kp_gpsw_oe.cpp $(RELIC_LIB_BLS_12_381) -lgmp
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
	$(CXX) -o objects/klin_kp_10000 $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) kLin_KP.cpp $(RELIC_LIB_BLS_12_381) -lgmp
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
	$(CXX) -o objects/abe_circuit_oe $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) abe_circuit_w_fanout_oe.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/abe_circuit_ok $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) abe_circuit_w_fanout_ok.cpp $(RELIC_LIB_BLS_12_381) -lgmp	
	$(CXX) -o objects/abe_circuit_base $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) abe_circuit_w_fanout_base.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/abe_circuit_gopt $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ) $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) abe_circuit_w_fanout_gopt.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/abe_circuit_pre $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ)  $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) abe_circuit_w_fanout_pre.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/abe_circuit_api $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ)  $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) abe_circuit_w_fanout_api.cpp $(RELIC_LIB_BLS_12_381) -lgmp
	$(CXX) -o objects/abe_circuit_GAP $(CXXFLAGS) $(ARITH_OBJ) $(POLICY_OBJ)  $(UTIL_OBJ) -I$(RELIC_INCLUDE_BLS_12_381) abe_circuit_w_fanout_GAP.cpp $(RELIC_LIB_BLS_12_381) -lgmp
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

