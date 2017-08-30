#ifndef SLRGEN_H
#define	SLRGEN_H

#include <iostream>
#include "test/testdata.h"
#include "index/PgSAIndexFactory.h"

using namespace PgSAReadsSet;
using namespace PgSAIndex;

namespace CLRgen {
	/**
     * Returne the neighbours of the k-mer kMer in the graph.
     * left: set to 1 to get the neighbours on the left of the k-mer,
     * and to 0 to get the neighbours of the right of the k-mer.
     */ 
    vector<StandardOccurrence> getNeighbours(string kMer, int left);
    
    /**
     * Extends the CLR LR on the left, on a maximum distance of extLen.
     */ 
    void extendLeft(int extLen, string &LR);
    
    /**
     * Extends the CLR LR on the right, on a maximum distance of extLen.
     */
    void extendRight(int extLen, string &LR);
    
    /**
     * Links src to dst by traversing the graph.
     * Parametrs:
     * src: source seed
     * dst: destination seed
     * curK: current k-mer size for this linking
     * visited: set of already visited seeds
     * curback: current number of backtracks
     * dist: current extension distance
     * curExt: current extension sequence
     * fRes: final result, the function updates it if src and dst could be linked
     */ 
    int link(string src, string dst, int curK, set<int> &visited, int* curback, int dist, string curExt, string &fRes);
    
    /**
	 * Links together the seeds contained in the file tolink, and outputs the corrected long read on the standard output.
	 * Parameters:
	 * index: PgSA index
	 * k: k-mer size
	 * tolink: files containing the seeds to be linked
	 * templateId: id of the original template long read
	 * seedsoverlap: minimum overlap required to allow the merging of two overlapping seeds
	 * backtracks: maximum number of backtracks allowed
	 * seedskips: maximum number of seeds that can be skipped
	 */ 
    void generateCLR(PgSAIndexStandard* index, int k, string tolink, string templateId, int seedsoverlap, int minoverlap, int backtracks, int seedskips);
}

#endif	/* SLRGEN_H */

