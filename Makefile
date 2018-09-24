CC=g++
CFLAGS  = -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11

PGSA_LIB=PgSA/dist/pgsalib/GNU-Linux-x86/
PGSA_SRC=PgSA/src/

LDFLAGS=-lPgSA -L$(PGSA_LIB) -lpthread

KMC_API_DIR = KMC/kmc_api
KMC_API_OBJS = \
$(KMC_API_DIR)/mmer.o \
$(KMC_API_DIR)/kmc_file.o \
$(KMC_API_DIR)/kmer_api.o
KMC_QUERY_DIR = src/kmc_query
KMC_QUERY_OBJS = \
$(KMC_QUERY_DIR)/nc_utils.o \
$(KMC_QUERY_DIR)/kmc_query.o

EXECS=CLRgen
all: $(EXECS)

CLRgen: seedsMerging.o seedsLinking.o $(KMC_QUERY_OBJS) $(KMC_API_OBJS) CLRgen.o
	$(CC) -o bin/CLRgen src/seedsMerging.o src/seedsLinking.o $(KMC_QUERY_OBJS) $(KMC_API_OBJS) src/CLRgen.o $(LDFLAGS) -Wl,-R$(PGSA_LIB)

seedsMerging.o: src/seedsMerging.cpp
	$(CC) -o src/seedsMerging.o -c src/seedsMerging.cpp $(CFLAGS)

seedsLinking.o: src/seedsLinking.cpp src/seedsMerging.h src/seedsLinking.h
	$(CC) -o src/seedsLinking.o -c src/seedsLinking.cpp $(CFLAGS) $(LDFLAGS) -I$(PGSA_SRC)

$(KMC_QUERY_OBJS): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

CLRgen.o: src/CLRgen.cpp src/seedsLinking.h
	$(CC) -o src/CLRgen.o -c src/CLRgen.cpp $(CFLAGS) $(LDFLAGS) -I$(PGSA_SRC)

clean:
	rm -Rf src/*.o src/kmc_query/*.o bin/CLRgen
