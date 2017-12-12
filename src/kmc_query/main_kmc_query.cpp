#include "stdafx.h"
#include <iostream>
#include "../kmc_api/kmc_file.h"
#include "nc_utils.h"
#include <iostream>
#include <chrono>

//~ using namespace std;

int _tmain(int argc, char* argv[]) {
	CKMCFile kmer_database;
	std::string file = string(argv[1]);
	
	kmer_database.OpenForRA(file);
	kmer_database.SetMinCount(1);
	kmer_database.SetMaxCount(10000000);
	
	uint32 _kmer_length;
	uint32 _mode;
	uint32 _counter_size;
	uint32 _lut_prefix_length;
	uint32 _signature_len;
	uint32 _min_count;
	uint64 _max_count;
	uint64 _total_kmers;
	kmer_database.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
	
	CKmerAPI kmer_object(_kmer_length);
	
	std::string line;
	uint64 counter;
	time_t t;
	while (std::getline(cin, line)) {
		auto start = std::chrono::high_resolution_clock::now();
		kmer_object.from_string(line);
		if (kmer_database.CheckKmer(kmer_object, counter)) {
			auto finish = std::chrono::high_resolution_clock::now();
			std::cout << "time : " << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() << std::endl;
			std::cout << line << " : " << counter << std::endl;
		} else {
			auto finish = std::chrono::high_resolution_clock::now();
			std::cout << "time : " << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() << std::endl;
		}
	}
	
	kmer_database.Close();
	
	return EXIT_SUCCESS;
}
