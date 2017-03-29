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
    int repeat = 11;
    int testKmersNumber = 1;
    vector<unsigned short> kLengths = {45};
    bool sFlag = false; // scramble reads (for uncorrecly concatenated pair-ended data)
    bool fFlag = false; // filter TTTTTT.....TTTT reads (for compatibility with CGk tests)
    bool pFlag = false; // query by position
    string cacheFile;
    string kParam;
    size_t pos;
    size_t found;
    char* tolink;
    char* tpl;
    int seedsoverlap;
    int minoverlap;
    int backtracks;
    int seedskips;

    while ((opt = getopt(argc, argv, "sfpr:n:c:k:l:t:o:m:b:i:?")) != -1) {
        switch (opt) {
        case 'r':
            repeat = atoi(optarg);
            break;
        case 'n':
            testKmersNumber = atoi(optarg);
            break;
        case 'c':
            cacheFile = optarg;
            break;
        case 'k':
            kLengths.clear();
            kParam = optarg;
            pos = 0;
            found = 0;
            while((found = kParam.find_first_of(",;|", pos)) != string::npos) {
                kLengths.push_back(atoi(kParam.substr(pos, found - pos).c_str()));
                pos = found + 1;
            }
            kLengths.push_back(atoi(kParam.substr(pos).c_str()));
            break;
        case 's':
            sFlag = true;
            break;
        case 'f':
            fFlag = true;
            break;
        case 'p':
            pFlag = true;
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
	case 'i':
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

    //cout << "*****************************************************************************\n";
    //cout << idx->getDescription();

    for(unsigned int i = 0; i < kLengths.size(); i++) {
        PgSATest::runTest(idx, repeat, testKmersNumber, kLengths[i], sFlag, fFlag, pFlag, kLengths[i], tolink, tpl, seedsoverlap, minoverlap, backtracks, seedskips);
    }
    delete(idx);

    exit(EXIT_SUCCESS);
}
