#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <algorithm> 
#include <string>

using namespace std;
 
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
	string seq;
	
	bool operator<(const seed_t& s2) const {
	  if (pos < s2.pos) {
		  return true;
	  } else if (pos == s2.pos && alen < s2.alen) {
		  return true;
	  } else if (pos == s2.pos && alen == s2.alen && matches < s2.matches) {
		  return true;
	  } else if (pos == s2.pos && alen == s2.alen && matches == s2.matches && seq < s2.seq) {
		  return true;
	  } else {
		  return false;
	  }
	}
};

/**
 * Computes the backtrack table of the string s. Auxiliary function used
 * to compute the overlap length between two strings.
 */
int* computeBacktrackTable(string s);

/**
 * Computes the overlap length between strings s1 and s2.
 */
int overlapLength(string s1, string s2);
 
/**
 * Merges the seeds contained in the vector seeds, if their alignment positions
 * indicate that they overlap over a greater length than minOverlap, and if 
 * their overlapping sequences match.
 */ 
void mergeOverlappingPosSeeds(vector<seed_t> &seeds, unsigned minOverlap);
 
/**
* Merges the seeds contained in the vector seeds, if their sequences overlap
* on a greater length than minOverlap.
*/ 
void mergeOverlappingSeqSeeds(vector<seed_t> &seeds, int minOverlap);

/**
 * Returns the length of the sequence of the template of identifier tpl.
 * len: length of the identifier.
 */ 
int getTemplateLength(string tpl);

/**
 * Reads the alignments stored in the file alFile, and 
 * returns a vector containing the corresponding seeds.
 */
vector<seed_t> readAlignmentFile(string alFile);
 
/**
 * Reads the alignments stored in the file alFile, and returns a vector
 * of the seeds, after the two merging steps, allowing a minimum overlap
 * of length minOverlap.
 */ 
vector<seed_t> processSeeds(string alFile, unsigned minOverlap);


