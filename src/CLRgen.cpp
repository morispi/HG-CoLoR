#include "index/cache/persistence/CountQueriesCachePersistence.h"
#include "suffixarray/persistence/SuffixArrayPersistence.h"
#include "seedsLinking.h"
#include "helper.h"
#include <stdlib.h>
#include <unistd.h>

PgSAIndexStandard* prepareIndex(string idxFile, string cacheFile) {
    return PgSAIndexFactory::getPgSAIndexStandard(idxFile, cacheFile, false);
}

int main(int argc, char *argv[]) {

    int opt; // current option
    int maxorder = 100;
    string cacheFile;
    string tmpDir;
    int seedsoverlap = maxorder - 1;
    int seedsdistance = 10;
    int minorder = maxorder / 2;
    int maxbranches = 1500;
    int seedskips = 5;
    int nbThreads = 1;
    int mismatches = 3;

    while ((opt = getopt(argc, argv, "c:K:t:o:d:k:b:s:j:m:?")) != -1) {
        switch (opt) {
			case 'c':
				cacheFile = optarg;
				break;
			case 'K':
				maxorder = atoi(optarg);
				break;
			case 't':
				tmpDir = optarg;
				break;
			case 'o':
				seedsoverlap = atoi(optarg);
				break;
			case 'd':
				seedsdistance = atoi(optarg);
				break;
			case 'k':
				minorder = atoi(optarg);
				break;
			case 'b':
				maxbranches = atoi(optarg);
				break;
			case 's':
				seedskips = atoi(optarg);
				break;
			case 'm':
				mismatches = atoi(optarg);
				break;
			case 'j':
				nbThreads = atoi(optarg);
				break;
			case '?':
			default: /* '?' */
				fprintf(stderr, "Usage: %s [-K maxOrder] [-t tmpDir] [-o seedsOverlap] [-d seedsMaxDistances] [-k minOrder] [-b maxBranches] [-s seedsSkips] [-j threadsNb] [-c cachefile] indexfile\n\n", argv[0]);
				exit(EXIT_FAILURE);
        }
    }

    if (optind != (argc - 1)) {
        fprintf(stderr, "%s: Expected only index name after options\n", argv[0]);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    string idxFile(argv[optind++]);

    PgSAIndexStandard* idx = prepareIndex(idxFile, cacheFile);
    
    CLRgen::startCorrection(idx, maxorder, tmpDir, seedsdistance, seedsoverlap, minorder, maxbranches, seedskips, mismatches, nbThreads);
    
    delete(idx);

    exit(EXIT_SUCCESS);
}
