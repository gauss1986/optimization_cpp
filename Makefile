CXXFLAGS =  -Wall -O2 -I/opt/intel/composer_xe_2013.3.171/mkl/include -I${SCINET_BOOST_INC} -I${SCINET_ARMADILLO_INC} -I${SCINET_R_BASE}/lib64/R/include/
LDFLAGS =  -L${MKLROOT}/lib/intel64 -L${SCINET_BOOST_LIB} -L${SCINET_R_LIB}/R/lib -L${SCINET_PCRE_LIB} -L${SCINET_XZ_LIB} -lR -llzma -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lgsl -lgslcblas -lpcre 
LDLIBS =  -lpthread -lboost_thread -lm
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

testR_main.o: testR_main.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o testR_main.o testR_main.cpp

testR.o: testR.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o testR.o testR.cpp

testR: testR.o testR_main.o
	$(CXX) $(LDFLAGS) -o testR testR.o testR_main.o $(LDLIBS)

clean:
	rm -f test maxsharpe maxsharpe.o txtIO.o stat.o testR.o testR_main.o testR
