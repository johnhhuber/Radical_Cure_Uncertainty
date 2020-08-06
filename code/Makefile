# define Makefile variables 
CC = g++
CFLAGS = -std=c++11 -O4

main.out: Source_Trial.o com.o Individual.o Intervention.o linpack.o Mosquito.o Params.o Population.o randlib.o Simulation.o Trial.o
	$(CC) $(CFLAGS) Source_Trial.o com.o Individual.o Intervention.o linpack.o Mosquito.o Params.o Population.o randlib.o Simulation.o Trial.o -o main.out
  
Source_Trial.o: Source_Trial.cpp Simulation.h randlib.h
	$(CC) $(CFLAGS) -c Source_Trial.cpp

Simulation.o: Simulation.cpp Simulation.h Intervention.h Trial.h
	$(CC) $(CFLAGS) -c Simulation.cpp

Trial.o: Trial.cpp Trial.h Individual.h Population.h Params.h
	$(CC) $(CFLAGS) -c Trial.cpp

Mosquito.o: Mosquito.cpp Population.h
	$(CC) $(CFLAGS) -c Mosquito.cpp

Intervention.o: Intervention.cpp Intervention.h Population.h
	$(CC) $(CFLAGS) -c Intervention.cpp

Population.o: Population.cpp Population.h Individual.h
	$(CC) $(CFLAGS) -c Population.cpp

Individual.o: Individual.cpp Individual.h Params.h
	$(CC) $(CFLAGS) -c Individual.cpp

Params.o: Params.cpp Params.h
	$(CC) $(CFLAGS) -c Params.cpp

com.o: com.cpp randlib.h
	$(CC) $(CFLAGS) -c com.cpp

randlib.o: randlib.cpp randlib.h
	$(CC) $(CFLAGS) -c randlib.cpp

linpack.o: linpack.cpp
	$(CC) $(CFLAGS) -c linpack.cpp

clean: 
	rm -f main.out com.o Individual.o Intervention.o linpack.o Mosquito.o Params.o Population.o randlib.o Simulation.o Trial.o 

 
