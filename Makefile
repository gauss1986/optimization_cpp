CXXFLAGS =  -Wall -O2 -I/opt/intel/composer_xe_2013.3.171/mkl/include -I${SCINET_BOOST_INC} -I${SCINET_ARMADILLO_INC} 
LDFLAGS =  -L${MKLROOT}/lib/intel64 -L${SCINET_BOOST_LIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lgsl -lgslcblas
LDLIBS = -lpthread -lboost_thread -lm
CXX = icpc

all: maxsharpe

txtIO.o: txtIO.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o txtIO.o txtIO.cpp

stat.o: stat.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o stat.o stat.cpp

maxsharpe.o: maxsharpe.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o maxsharpe.o maxsharpe.cpp

maxsharpe: txtIO.o stat.o maxsharpe.o
	$(CXX) $(LDFLAGS) -o maxsharpe txtIO.o stat.o maxsharpe.o $(LDLIBS)

test: test.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) 

clean:
	rm -f test maxsharpe maxsharpe.o txtIO.o stat.o
