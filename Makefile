CFLAGS = -I/usr/include/eigen3 
all: main.cpp
	g++ $(CFLAGS) -Wall -O2 -o smallpt4d main.cpp -pthread -std=c++11 -fopenmp