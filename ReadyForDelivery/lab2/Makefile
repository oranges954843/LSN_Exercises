CFLAGS = -Wall -O3 --std=c++11


main1.exe: main1.o integrale.o random.o integranda.o
	g++ main1.o integrale.o random.o integranda.o -o main1.exe

main2.exe: main2.o random.o 
	g++ main2.o random.o -o main2.exe 

main1.o: main1.cpp integrale.h random.h integranda.h FunzioneBase.h  
	g++ -c main1.cpp -o main1.o $(CFLAGS)

main2.o: main2.cpp random.h  
	g++ -c main2.cpp -o main2.o $(CFLAGS)

integrale.o: integrale.cpp integranda.h integrale.h FunzioneBase.h random.h
	g++ -c integrale.cpp -o integrale.o $(CFLAGS)

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)

integranda.o: integranda.cpp integranda.h FunzioneBase.h
	g++ -c integranda.cpp -o integranda.o 

exe1:
	./main1.exe

exe2:
	./main2.exe

clean:
	rm *.o *.exe seed.out

clean~:
	rm *~
