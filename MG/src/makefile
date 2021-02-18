CC = g++
IDIR =/usr/local/include
LDIR =/usr/local/lib
current_dir = $(shell pwd)
CFLAGS = -Wall -ansi -O3 -DCPU_ONLY -std=c++11 -I$(IDIR) -I$(current_dir)/../include/ -c
LFLAGS = -L$(LDIR) -lm  -lgsl

all: evolution_degeneracy, lane_emden

evolution_degeneracy: evolution_degeneracy.o
	$(CC) $(LFLAGS) -o evolution_degeneracy evolution_degeneracy.o

evolution_degeneracy.o: evolution_degeneracy.cpp
	$(CC) $(CFLAGS) evolution_degeneracy.cpp

lane_emden: lane_emden.o
	$(CC) $(LFLAGS) -o lane_emden lane_emden.o

lane_emden.o: lane_emden.cpp
	$(CC) $(CFLAGS) lane_emden.cpp

clean:
	rm *.o evolution_degeneracy, lane_emden
