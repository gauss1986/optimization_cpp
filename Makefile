CXXFLAGS =  -O1 -I/opt/intel/composer_xe_2013.3.171/mkl/include -I${SCINET_BOOST_INC} -I${SCINET_ARMADILLO_INC}
LDFLAGS =  -L${MKLROOT}/lib/intel64 -L${SCINET_BOOST_LIB} -L${SCINET_ARMADILLO_LIB} -larmadillo -llzma -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lgsl -lgslcblas 
LDLIBS =  -lpthread -lboost_thread -lm
CXX = icpc

all: main

ticktock.o: ticktock.cc
	$(CXX) -c $(CXXFLAGS) -I . -o ticktock.o ticktock.cc

txtIO.o: txtIO.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o txtIO.o txtIO.cpp

stat.o: stat.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o stat.o stat.cpp

OLS.o: OLS.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o OLS.o OLS.cpp

maxshp.o: maxshp.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o maxshp.o maxshp.cpp

comp_shp.o: comp_shp.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o comp_shp.o comp_shp.cpp

main.o: main.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o main.o main.cpp

main: ticktock.o txtIO.o stat.o OLS.o maxshp.o comp_shp.o main.o
	$(CXX) $(LDFLAGS) -o main ticktock.o txtIO.o stat.o OLS.o maxshp.o comp_shp.o main.o $(LDLIBS)

test: test.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) 

clean:
	rm -f test main ticktock.o main.o txtIO.o stat.o OLS.o maxshp.o comp_shp.o
