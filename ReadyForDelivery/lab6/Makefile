CC = g++
CFLAGS = -Wall -O3 --std=c++11

Metropolis_ISING_1D.exe : Metropolis_ISING_1D.o random.o
	$(CC) random.o Metropolis_ISING_1D.o -o Metropolis_ISING_1D.exe
Metropolis_ISING_1D.o : Metropolis_ISING_1D.cpp
	$(CC) -c Metropolis_ISING_1D.cpp -o Metropolis_ISING_1D.o $(CFLAGS)
Gibbs_ISING_1D.exe : Gibbs_ISING_1D.o random.o
	$(CC) random.o Gibbs_ISING_1D.o -o Gibbs_ISING_1D.exe
Gibbs_ISING_1D.o : Gibbs_ISING_1D.cpp
	$(CC) -c Gibbs_ISING_1D.cpp -o Gibbs_ISING_1D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o *.out *.exe

clean_metro :
	rm *.o *.metro *.out *.exe

clean_gibbs :
	rm *.o *.gibbs *.out *.exe
