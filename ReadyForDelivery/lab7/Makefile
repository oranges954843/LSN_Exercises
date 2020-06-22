CC = g++
CFLAGS = -Wall -O3 --std=c++11

Istant_Values.exe : Istant_Values.o random.o
	$(CC) random.o Istant_Values.o -o Istant_Values.exe
Istant_Values.o : Istant_Values.cpp Monte_Carlo_NVT.h
	$(CC) -c Istant_Values.cpp -o Istant_Values.o $(CFLAGS)

Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o random.o
	$(CC) random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp Monte_Carlo_NVT.h
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.0 *.o *.exe *.out *.final

clean_ist :
	rm *.istant *.o *.exe *.out *.final

MolDyn_NVE.exe: MolDyn_NVE.o 
	g++ MolDyn_NVE.o -o MolDyn_NVE.exe

MolDyn_NVE.o: MolDyn_NVE.cpp MolDyn_NVE.h
	g++ -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)

clean_NVE:
	rm *.o *.out *.exe *.final *.istant
