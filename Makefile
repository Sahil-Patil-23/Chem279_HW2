# GNU C++ Compiler
CPP  =  g++

CPP_FLAGS  =  -std=c++11 -o

SRC_DIR  =  src
BIN_DIR  = bin


Numerical_Integration:
	cd $(SRC_DIR) && $(CPP) $(CPP_FLAGS) ../$(BIN_DIR)/HW2_Q1 HW2_Numerical.cpp

Run_Numerical_Integration:
	cd $(BIN_DIR) && ./HW2_Q1

Analytical_Integration:
	cd $(SRC_DIR) && $(CPP) $(CPP_FLAGS) ../$(BIN_DIR)/HW2_Q2 HW2_Analytical.cpp

Run_Analytical_Integration:
	cd $(BIN_DIR) && ./HW2_Q2

P1:
	make Numerical_Integration
	make Run_Numerical_Integration

P2:
	make Analytical_Integration
	make Run_Analytical_Integration

all:
	make P1
	make P2