#include "index/cache/persistence/CountQueriesCachePersistence.h"
#include "suffixarray/persistence/SuffixArrayPersistence.h"
#include "seedsLinking.h"

#include "helper.h"
#include <stdlib.h>    /* for exit */
#include <unistd.h>

PgSAIndexStandard* prepareIndex(string idxFile, string cacheFile) {

    return PgSAIndexFactory::getPgSAIndexStandard(idxFile, cacheFile, false);

};

int main(int argc, char *argv[])
{

    int opt; // current option
    int k = 64;
    string cacheFile;
    string kParam;
    size_t pos;
    size_t found;
    string tolink;
    string tpl;
    int seedsoverlap;
    int minoverlap;
    int backtracks;
    int seedskips;

    while ((opt = getopt(argc, argv, "c:k:l:t:o:m:b:s:?")) != -1) {
        switch (opt) {
			case 'c':
				cacheFile = optarg;
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'l':
				tolink = optarg;
				break;
			case 't':
				tpl = optarg;
				break;
			case 'o':
				seedsoverlap = atoi(optarg);
				break;
			case 'm':
				minoverlap = atoi(optarg);
				break;
			case 'b':
				backtracks = atoi(optarg);
				break;
			case 's':
				seedskips = atoi(optarg);
				break;
			case '?':
			default: /* '?' */
				fprintf(stderr, "Usage: %s [-k length] [-r no of repeats] [-n no of testkmers] [-c cachefile] [-p] [-s] [-f] indexfile\n\n",
						argv[0]);
				fprintf(stderr, "-p query by position\n-s scramble reads (for uncorrecly concatenated pair-ended data)\n-f -filter TTTTTT.....TTTT reads (for compatibility with CGk tests)\n\n");
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
    
    CLRgen::generateCLR(idx, k, tolink, tpl, seedsoverlap, minoverlap, backtracks, seedskips);
    
    delete(idx);

    exit(EXIT_SUCCESS);
}
