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

EXECS=HG-CoLoR
all: $(EXECS)

HG-CoLoR: seedsMerging.o seedsLinking.o reverseComplement.o $(KMC_QUERY_OBJS) $(KMC_API_OBJS) main.o
	$(CC) -o bin/HG-CoLoR src/reverseComplement.o src/seedsMerging.o src/seedsLinking.o $(KMC_QUERY_OBJS) $(KMC_API_OBJS) src/main.o $(LDFLAGS) -Wl,-R$(PGSA_LIB)

reverseComplement.o: src/reverseComplement.cpp
	$(CC) -o src/reverseComplement.o -c src/reverseComplement.cpp $(CFLAGS)

seedsMerging.o: src/seedsMerging.cpp
	$(CC) -o src/seedsMerging.o -c src/seedsMerging.cpp $(CFLAGS)

seedsLinking.o: src/seedsLinking.cpp src/seedsMerging.h src/seedsLinking.h src/reverseComplement.h
	$(CC) -o src/seedsLinking.o -c src/seedsLinking.cpp $(CFLAGS) $(LDFLAGS) -I$(PGSA_SRC)

$(KMC_QUERY_OBJS): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

main.o: src/main.cpp src/seedsLinking.h
	$(CC) -o src/main.o -c src/main.cpp $(CFLAGS) $(LDFLAGS) -I$(PGSA_SRC)

clean:
	rm -Rf src/*.o src/kmc_query/*.o bin/HG-CoLoR
