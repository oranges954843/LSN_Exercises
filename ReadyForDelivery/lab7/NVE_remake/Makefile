CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe: MolDyn_NVE.o 
	g++ MolDyn_NVE.o -o MolDyn_NVE.exe

MolDyn_NVE.o: MolDyn_NVE.cpp MolDyn_NVE.h
	g++ -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)

clean_NVE:
	rm *.o *.out *.exe *.final *.istant
