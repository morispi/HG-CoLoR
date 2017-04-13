#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <utility>
#include <vector>
#include <algorithm> 
#include <string>
 
/**
 * Data structure used to sore informations about the alignment of a seed
 * on a template.
 * alen: alignment length
 * pos: alignment position
 * tlen: template length
 * matches: number of matching bases
 * seq: DNA sequence of the alignment
 */ 
struct seed_t {
	int pos;
	int alen;
	int tlen;
	int matches;
	std::string seq;
	
	bool operator<(const seed_t& s2) const {
	  return pos < s2.pos;
	}
};

/**
 * Computes the backtrack table of the string s. Auxiliary function used
 * to compute the overlap length between two strings.
 */
int* computeBacktrackTable(std::string s);

/**
 * Computes the overlap length between strings s1 and s2.
 */
int overlapLength(std::string s1, std::string s2);
 
/**
 * Merges the seeds contained in the vector seeds, if their alignment positions
 * indicate that they overlap over a greater length than minOverlap, and if 
 * their overlapping sequences match.
 */ 
void mergeOverlappingPosSeeds(std::vector<seed_t> &seeds, unsigned minOverlap);
 
/**
* Merges the seeds contained in the vector seeds, if their sequences overlap
* on a greater length than minOverlap.
*/ 
void mergeOverlappingSeqSeeds(std::vector<seed_t> &seeds, int minOverlap);

/**
 * Returns the length of the sequence of the template of identifier tpl.
 * len: length of the identifier.
 */ 
int getTemplateLength(char* tpl, int len);

/**
 * Reads the alignments stored in the file alFile, and 
 * returns a vector containing the corresponding seeds.
 */
std::vector<seed_t> readAlignmentFile(char* alFile);
 
/**
 * Reads the alignments stored in the file alFile, and returns a vector
 * of the seeds, after the two merging steps, allowing a minimum overlap
 * of length minOverlap.
 */ 
std::vector<seed_t> processSeeds(char* alFile, unsigned minOverlap);


