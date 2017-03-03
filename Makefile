CXXFLAGS = -Wall -g -O3
LDFLAGS = -g
LDLIBS = 
CXX = g++

all: maxsharpe

maxsharpe.o: maxsharpe.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o maxsharpe.o maxsharpe.cpp

readtxt.o: readtxt.cpp
	$(CXX) -c $(CXXFLAGS) -I . -o readtxt.o readtxt.cpp

maxsharpe: maxsharpe.o readtxt.o
	$(CXX) $(LDFLAGS) -o maxsharpe maxsharpe.o readtxt.o $(LDLIBS)

clean:
	rm -f maxsharpe.o readtxt.o
