CXXFLAGS =  -O1 -I/opt/intel/composer_xe_2013.3.171/mkl/include -I${SCINET_BOOST_INC} -I${SCINET_ARMADILLO_INC}
LDFLAGS =  -L${MKLROOT}/lib/intel64 -L${SCINET_BOOST_LIB} -L${SCINET_PCRE_LIB} -L${SCINET_XZ_LIB} -L${SCINET_ARMADILLO_LIB} -larmadillo -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread 
LDLIBS =  -lpthread -lboost_thread -lm
CXX = icpc

all: main

ticktock.o: ticktock.cc
	$(CXX) -c $(CXXFLAGS) -I . -o ticktock.o ticktock.cc

misc.o: misc.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o misc.o misc.cpp

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

main_real.o: main_real.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o main_real.o main_real.cpp

main: ticktock.o misc.o txtIO.o stat.o OLS.o maxshp.o main.o
	$(CXX) $(LDFLAGS) -o main ticktock.o misc.o txtIO.o stat.o OLS.o maxshp.o main.o $(LDLIBS)

real: ticktock.o misc.o txtIO.o stat.o OLS.o maxshp.o main_real.o
	$(CXX) $(LDFLAGS) -o real ticktock.o misc.o txtIO.o stat.o OLS.o maxshp.o main_real.o $(LDLIBS)

clean:
	rm -f main real ticktock.o misc.o main.o main_real.o txtIO.o stat.o OLS.o maxshp.o 
