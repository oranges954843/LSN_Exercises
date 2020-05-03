CFLAGS = -Wall -O3 --std=c++11


main1.exe: main1.o random.o
	g++ main1.o random.o -o main1.exe 

main1.o: main1.cpp random.h  
	g++ -c main1.cpp -o main1.o $(CFLAGS)

main2.exe: main2.o random.o
	g++ main2.o random.o -o main2.exe  

main2.o: main2.cpp random.h  
	g++ -c main2.cpp -o main2.o $(CFLAGS)

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)

exe1:
	./main1.exe

exe2:
	./main2.exe

clean:
	rm *.o *.exe seed.out

