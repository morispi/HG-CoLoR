#include "pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.h"
#include "index/cache/generator/CountQueriesCacheGenerator.h"
#include "index/cache/persistence/CountQueriesCachePersistence.h"
#include "suffixarray/generator/SuffixArrayGenerator.h"
#include "suffixarray/persistence/SuffixArrayPersistence.h"

#include "helper.h"
#include <stdlib.h>    /* for exit */
#include <unistd.h>

PseudoGenomeBase* preparePg(string srcFile, string pairFile) {
    PseudoGenomeBase* pgb = 0;
    
    if (PseudoGenomePersistence::isValidPseudoGenome(srcFile))
        pgb = PseudoGenomePersistence::readPseudoGenome(srcFile);
    else if (SuffixArrayPersistence::isValidSuffixArray(srcFile)) 
        pgb = SuffixArrayPersistence::readPgOnly(srcFile);
    else {         
        PseudoGenomeGeneratorFactory* pggf = new GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory();
        PseudoGenomeGeneratorBase* pggb = pggf->getGenerator(srcFile, pairFile);
        pgb = pggb->generatePseudoGenomeBase();        
        delete(pggb); 
        delete(pggf);
        pgb->buildRepetitiveReadsFilter();
    }
    
    if (pgb == 0) {
        fprintf(stderr, "Failed generating Pg\n");
        exit(EXIT_FAILURE);    
    }
    
    return pgb;
}

CountQueriesCacheBase* generateCache(PseudoGenomeBase* pgb) {
    CountQueriesCacheBase* cqcb = 0;

    cqcb = CountQueriesCacheGenerator::generateCountQueriesCache(pgb);
    
    return cqcb;
};

SuffixArrayBase* generateSA(PseudoGenomeBase* pgb, int rate, int fixed_min_k) {
    SuffixArrayBase* sab;

    if (rate == 1) 
        sab = SuffixArrayGenerator::generateDefaultSuffixArray(pgb, fixed_min_k);
    else if (rate > 1)
        sab = SuffixArrayGenerator::generateSparseSuffixArray(pgb, rate, fixed_min_k);    
    else {
        fprintf(stderr, "Unsupported compression rate: %d\n", rate);
        exit(EXIT_FAILURE);
    }
    return sab;
};


int main(int argc, char *argv[])
{

    int opt; // current option
    int compressionRate = 1; 
    int fixed_min_k = 1;
    bool pFlag = false;
    bool cFlag = false;

    while ((opt = getopt(argc, argv, "r:k:pc?")) != -1) {
        switch (opt) {
        case 'p':
            pFlag = true;
            break;
        case 'c':
            cFlag = true;
            break;
        case 'r':
            compressionRate = atoi(optarg);
            break;
        case 'k':
            fixed_min_k = atoi(optarg);
            break;
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-r rate] [-k fixed_k] [-c] [-p] readssrcfile [pairsrcfile] indexprefix\n\n",
                    argv[0]);
            fprintf(stderr, "-r compression rate [1 - 6]\n-k fixed kmer length \n-c generate cache file\n-p only Pg, no SA\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 2) || optind < (argc - 3)) {
        fprintf(stderr, "%s: Expected 2 or 3 arguments after options (found %d)\n", argv[0], argc-optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }
        
    string srcFile(argv[optind++]);
    string pairFile = "";
    if (optind == argc - 2)
        pairFile = argv[optind++];
    string idxPrefix(argv[optind++]);
    
    PseudoGenomeBase* pgb = preparePg(srcFile, pairFile);
    
    if (cFlag) {
        CountQueriesCacheBase* cqcb = generateCache(pgb);
        if (cqcb == 0) {
            fprintf(stderr, "Error: cache not available.\n");
        } else {
            CountQueriesCachePersistence::writeCache(cqcb, idxPrefix); 
            delete(cqcb);
        }
    }
        
    if (pFlag) {        
        PseudoGenomePersistence::writePseudoGenome(pgb, idxPrefix);
        delete(pgb);
    } 
    
    if (!pFlag && !cFlag) {
        SuffixArrayBase* sab = generateSA(pgb, compressionRate, fixed_min_k);
        SuffixArrayPersistence::writePgSA(sab, idxPrefix);
        delete(sab);
    }
    
//    cout << "DONE \n";
//    cin.ignore(1);    
    
    exit(EXIT_SUCCESS);
}
