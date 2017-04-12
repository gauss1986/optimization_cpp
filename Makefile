CXXFLAGS =  -O1 -g -ggdb -pg -I/usr/include/boost/ -I/usr/local/include -I/opt/intel/mkl/include  
LDLIBS =  -L/usr/local/lib64/ -L/usr/lib64/ -L/opt/intel/mkl/lib/intel64 -L/usr/lib64/boost 
LDFLAGS =  -lpthread -lboost_thread-mt -lm -lmkl_rt -llzma -lpthread -lgsl -lgslcblas -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -larmadillo -pg
CXX = gcc

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
