//LIBS:=`root-config --libs`
//INCS:=`root-config --cflags`
CC = g++
CFLAGS = -Wall -O3 --std=c++11

%.o : %.cpp
	g++ -Wall -c $< ${INCS}

qmc1d: qmc1d.o
	g++ -O3 -Wall -o $@ $^ ${LIBS}

VAR_MC_1DP.exe: VAR_MC_1DP.o random.o
	$(CC) random.o VAR_MC_1DP.o -o VAR_MC_1DP.exe

VAR_MC_1DP.o: VAR_MC_1DP.cpp VAR_MC_1DP.h
	$(CC) -c VAR_MC_1DP.cpp -o VAR_MC_1DP.o $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean_qmc:
	rm *.o qmc1d potential.dat kinetic.dat probability.dat

clean:
	rm *.o *.exe

clean_results:
	rm *.o *.exe *.out *.final

