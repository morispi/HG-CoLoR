#ifndef SEEDSLNK_H
#define	SEEDSLNK_H

#include <iostream>
#include "test/testdata.h"
#include "index/PgSAIndexFactory.h"
#include "../CTPL/ctpl_stl.h"
#include "seedsMerging.h"

using namespace PgSAReadsSet;
using namespace PgSAIndex;

// namespace CLRgen {
	/**
     * Returns the neighbours of the k-mer kMer in the de Bruijn graph of order k.
     * left: set to 1 to get the neighbours on the left of the k-mer,
     * and to 0 to get the neighbours on the right of the k-mer.
     */ 
    vector<string> getNeighbours(string kMer, int left);
    
    /**
     * Extends the CLR LR on the left, on a maximum distance of extLen.
     * Returns: The actual size of the extension.
     */ 
    unsigned extendLeft(unsigned extLen, string &LR);
    
    /**
     * Extends the CLR LR on the right, on a maximum distance of extLen.
     * Returns: The actual size of the extension.
     */
    unsigned extendRight(unsigned extLen, string &LR);
    
    /**
     * Attempts to link srcSeed to tgtSeed by traversing the graph.
     * Parameters:
     * srcSeed: source seed
     * tgtSeed: target seed
     * curK: order of the current de Bruijn graph
     * visited: set of already visited nodes
     * curBranches: current number of branches explorations
     * dist: current extension distance
     * curExt: current extension sequence
     * missingPart: string corresponding to the missing part of the long read. Updated if srcSeed and tgtSeed can be linked.
     * LRLen: length of the original long read
     * Returns: 1 if srcSeed can be linked to tgtSeed, 0 otherwhise.
     */ 
     int link(string srcSeed, string tgtSeed, unsigned curK, set<string> &visited, unsigned* curBranches, unsigned dist, string curExt, string &missingPart, unsigned LRLen);
    
    /**
	 * Generates the corrected long read for tpl from the seeds vector.
	 */ 
    std::pair<string, string> correctRead(int id, string tpl, vector<seed_t> seeds);
    
    /**
     * Launches the correction prodecure.
     * Parameters:
	 * index: PgSA index
	 * maxOrder: maximum order of the variable-order de Bruijn graph
	 * tmpDir: directory containing the temporary files (long readss to be corrected and associated seeds)
	 * seedsdistance: maximum distance allowed between two consecutive seeds during the second merging phase
	 * seedsoverlap: minimum overlap required to allow the merging of two overlapping seeds
	 * minorder: minimum order of the variable-order de Bruijn graph
	 * maxbranches: maximum number of branches exploration allowed
	 * maxseedsskips: maximum number of seeds that can be skipped
	 * mismatches: mismatch threshold tolerance
	 * nbThreads: number of threads to use
	 * longReadsFile : file containing the raw long reads to correct
	 * alignmentsFile : file containing the alignments of the SRs to the LRs
	 */ 
    void startCorrection(PgSAIndexStandard* index, unsigned maxorder, string tmpdir, unsigned seedsdistance, unsigned seedsoverlap, unsigned minorder, unsigned maxbranches, unsigned maxseedsskips, unsigned mismatches, unsigned nbThreads, string longReadsFile, string alignmentsFile);
// }

#endif	/* SEEDSLNK_H */
