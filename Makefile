CXXFLAGS =  -O1 -I/opt/intel/composer_xe_2013.3.171/mkl/include -I${SCINET_BOOST_INC} -I${SCINET_ARMADILLO_INC} -I${SCINET_R_BASE}/lib64/R/include/
LDFLAGS =  -L${MKLROOT}/lib/intel64 -L${SCINET_BOOST_LIB} -L${SCINET_R_LIB}/R/lib -L${SCINET_PCRE_LIB} -L${SCINET_XZ_LIB} -lR -llzma -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lgsl -lgslcblas -lpcre 
LDLIBS =  -lpthread -lboost_thread -lm
CXX = icpc

all: main

txtIO.o: txtIO.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o txtIO.o txtIO.cpp

stat.o: stat.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o stat.o stat.cpp

OLS.o: OLS.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o OLS.o OLS.cpp

maxshp.o: maxshp.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o maxshp.o maxshp.cpp

main.o: main.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o main.o main.cpp

main: txtIO.o stat.o OLS.o maxshp.o main.o
	$(CXX) $(LDFLAGS) -o main txtIO.o stat.o OLS.o maxshp.o main.o $(LDLIBS)

test: test.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) 

testR_main.o: testR_main.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o testR_main.o testR_main.cpp

testR.o: testR.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o testR.o testR.cpp

testR: testR.o testR_main.o
	$(CXX) $(LDFLAGS) -o testR testR.o testR_main.o $(LDLIBS)

clean:
	rm -f test main main.o txtIO.o stat.o testR.o testR_main.o OLS.o maxshp.o testR
