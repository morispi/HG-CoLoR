CC=g++
CFLAGS=-std=c++11 -O3
PGSA_LIB=$(PGSA_PATH)dist/pgsalib/GNU-Linux-x86/
PGSA_SRC=$(PGSA_PATH)src
LDFLAGS=-lPgSA -L$(PGSA_LIB)
EXECS=PgSAgen CLRgen
all: $(EXECS)

CLRgen: seedsMerging.o seedsLinking.o CLRgen.o
	$(CC) -o bin/CLRgen src/seedsMerging.o src/seedsLinking.o src/CLRgen.o $(LDFLAGS) -Wl,-R$(PGSA_LIB)

seedsMerging.o: src/seedsMerging.cpp
	$(CC) -o src/seedsMerging.o -c src/seedsMerging.cpp $(CFLAGS)

seedsLinking.o: src/seedsLinking.cpp src/seedsMerging.h src/seedsLinking.h
	$(CC) -o src/seedsLinking.o -c src/seedsLinking.cpp $(CFLAGS) $(LDFLAGS) -I$(PGSA_SRC)

CLRgen.o: src/CLRgen.cpp src/seedsLinking.h
	$(CC) -o src/CLRgen.o -c src/CLRgen.cpp $(CFLAGS) $(LDFLAGS) -I$(PGSA_SRC)

PgSAgen: PgSAgen.o
	$(CC) -o bin/PgSAgen src/PgSAgen.o $(LDFLAGS) -Wl,-R$(PGSA_LIB)

PgSAgen.o: src/PgSAgen.cpp
	$(CC) -o src/PgSAgen.o -c src/PgSAgen.cpp $(CFLAGS) $(LDFLAGS) -I$(PGSA_SRC)

clean:
	rm -rf src/*.o
