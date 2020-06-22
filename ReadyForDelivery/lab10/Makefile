CC = g++
CFLAGS = -Wall -O3 --std=c++11


SA_TSP.exe : SA_TSP.o random.o
	$(CC) SA_TSP.o random.o -o SA_TSP.exe

SA_TSP.o : SA_TSP.cpp 
	$(CC) -c SA_TSP.cpp -o SA_TSP.o $(CFLAGS) 

GA_TSP.exe : GA_TSP.o random.o
	$(CC) GA_TSP.o random.o -o GA_TSP.exe

GA_TSP.o : GA_TSP.cpp
	$(CC) -c GA_TSP.cpp -o GA_TSP.o $(CFLAGS) 

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o *.out *.exe

