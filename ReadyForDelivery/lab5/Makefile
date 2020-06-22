CFLAGS = -Wall -O3 --std=c++11


main.exe: main.o random.o
	g++ main.o random.o -o main.exe 

main.o: main.cpp random.h main.h 
	g++ -c main.cpp -o main.o $(CFLAGS)

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)

exe:
	./main.exe

clean:
	rm *.o *.exe *.out *.dat

