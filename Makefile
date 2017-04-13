CC=g++
CFLAGS=-std=c++11 -O3
PGSA_LIB=$(PGSA_PATH)dist/pgsalib/GNU-Linux-x86/
PGSA_SRC=$(PGSA_PATH)src
LDFLAGS=-lPgSA -L$(PGSA_LIB)
EXECS=PgSAgen SLRgen
all: $(EXECS)

SLRgen: seedsMerging.o seedsLinking.o SLRgen.o
	$(CC) -o SLRgen src/seedsMerging.o src/seedsLinking.o src/SLRgen.o $(LDFLAGS) -Wl,-R$(PGSA_LIB)

seedsMerging.o: src/seedsMerging.cpp
	$(CC) -o src/seedsMerging.o -c src/seedsMerging.cpp $(CFLAGS)

seedsLinking.o: src/seedsLinking.cpp src/seedsMerging.h src/seedsLinking.h
	$(CC) -o src/seedsLinking.o -c src/seedsLinking.cpp $(CFLAGS) $(LDFLAGS) -I$(PGSA_SRC)

SLRgen.o: src/SLRgen.cpp src/seedsLinking.h
	$(CC) -o src/SLRgen.o -c src/SLRgen.cpp $(CFLAGS) $(LDFLAGS) -I$(PGSA_SRC)

PgSAgen: PgSAgen.o
	$(CC) -o PgSAgen src/PgSAgen.o $(LDFLAGS) -Wl,-R$(PGSA_LIB)

PgSAgen.o: src/PgSAgen.cpp
	$(CC) -o src/PgSAgen.o -c src/PgSAgen.cpp $(CFLAGS) $(LDFLAGS) -I$(PGSA_SRC)

clean:
	rm -rf src/*.o
