CC = g++
IDIR =/usr/local/include
LDIR =/usr/local/lib
current_dir = $(shell pwd)
CFLAGS = -Wall -ansi -O3 -DCPU_ONLY -std=c++11 -I$(IDIR) -I$(current_dir)/../include/ -c
LFLAGS = -L$(LDIR) -lm  -lgsl

all: evolution_degeneracy_MMSM

evolution_degeneracy_MMSM: evolution_degeneracy_MMSM.o
	$(CC) $(LFLAGS) -o evolution_degeneracy_MMSM evolution_degeneracy_MMSM.o

evolution_degeneracy_MMSM.o: evolution_degeneracy_MMSM.cpp
	$(CC) $(CFLAGS) evolution_degeneracy_MMSM.cpp

clean:
	rm *.o evolution_degeneracy evolution_degeneracy_MMSM
