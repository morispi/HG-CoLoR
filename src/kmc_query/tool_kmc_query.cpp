#include "stdafx.h"
#include <iostream>
#include "../kmc_api/kmc_file.h"
#include "nc_utils.h"

CKMCFile database;
CKmerAPI kmer_object;


void openDatabase(std::string file) {
	if (!database.OpenForRA(file)) {
		std::cerr << "Can't open k-mer database file \"file\"." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	CKmerAPI k(database.KmerLength());
	kmer_object = k;
}

void closeDatabase() {
	if (!database.Close()) {
		std::cerr << "Error closing k-mer database file \"file\"." << std::endl;
		exit(EXIT_FAILURE);
	}
}

uint64 getOccNb(string kmer) {
	uint64 counter;
	
	kmer_object.from_string(kmer);
	if (!database.CheckKmer(kmer_object, counter)) {
		counter = 0;
	}
	
	return counter;
}
