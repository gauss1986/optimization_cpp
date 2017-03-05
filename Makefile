CXXFLAGS =  -Wall -O3 -I/opt/intel/composer_xe_2013.3.171/mkl/include 
LDFLAGS =  -L/opt/intel/composer_xe_2013.3.171/mkl/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lgsl -lgslcblas 
LDLIBS = -lpthread -lm
CXX = icc

all: maxsharpe

maxsharpe.o: maxsharpe.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o maxsharpe.o maxsharpe.cpp

txtIO.o: txtIO.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o txtIO.o txtIO.cpp

maxsharpe: maxsharpe.o txtIO.o
	$(CXX) $(LDFLAGS) -o maxsharpe maxsharpe.o txtIO.o $(LDLIBS)

clean:
	rm -f maxsharpe maxsharpe.o txtIO.o
