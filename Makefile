## module load intel/13.1.1
# module load armadillo/3.910.0
# # #
# #
CXX=icpc
# #
EXTRA_LIB_FLAGS =  -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
# #
LIB_FLAGS = $(EXTRA_LIB_FLAGS)
# #
OPT = -O3
# # ## As the Armadillo library uses recursive templates,
# # ## compilation times depend on the level of optimisation:
# # ##
# # ## -O0: quick compilation, but the resulting program will be slow
# # ## -O1: good trade-off between compilation time and execution speed
# # ## -O2: produces programs which have almost all possible speedups,
# # ##      but compilation takes longer
# #
# #
# #
CXXFLAGS = $(OPT) -I$(SCINET_ARMADILLO_INC)
# #
all: armadillotest
# #
# #
armadillotest: test.cpp
	$(CXX) $(CXXFLAGS) $(EXTRAFLAGS)  -o $@  $<  $(LIB_FLAGS)
# 	#
# 	#
.PHONY: clean
# 	#
clean:
# 	#
# 	#

