CXXFLAGS =  -Wall -O3 -I/opt/intel/composer_xe_2013.3.171/mkl/include -I${SCINET_BOOST_INC} 
LDFLAGS =  -L/opt/intel/composer_xe_2013.3.171/mkl/lib -L${SCINET_BOOST_LIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lgsl -lgslcblas  
LDLIBS = -lpthread -lboost_thread -lm
CXX = icc

all: maxsharpe

txtIO.o: txtIO.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o txtIO.o txtIO.cpp

maxsharpe.o: maxsharpe.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o maxsharpe.o maxsharpe.cpp

maxsharpe: txtIO.o maxsharpe.o
	$(CXX) $(LDFLAGS) -o maxsharpe txtIO.o maxsharpe.o $(LDLIBS)

clean:
	rm -f maxsharpe maxsharpe.o txtIO.o
