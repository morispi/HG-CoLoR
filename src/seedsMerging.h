#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <utility>
 
/**
 * Data structure to sore informations about the alignment of a SR
 * on a template.
 * alen: alignment length
 * pos: alignment position
 * tlen: template length
 * matches: number of matching bases
 * seq: DNA sequence of the alignment
 */ 
typedef struct seed_s {
	int alen;
	int pos;
	int tlen;
	int matches;
	char* seq;
} seed_t; 

/**
 * Turns the string str of length len in uppercase.
 */ 
char* strToUpper(char* str, int len);

/**
 * Computes the backtrack table of the string s.
 */
int* computeBacktrackTable(char* s);

/**
 * Computes the overlap length between strings s1 and s2.
 */
int overlapLength(char* s1, char* s2);

/**
 * Compares the alignment positions of seeds a1 and a2.
 * Returns a negative value if a1 aigned on the left of a2,
 * zero if they aligned at the same position,
 * and a positive value if a1 aligned on the right of a2.
 */ 
int compareAlignmentsPos(const void* r1, const void* r2);

/**
 * Merges the seeds contained in the structure of size *nb, pointed
 * by seeds if their alignment positions indicate that they overlap
 * on a greater length than minOverap.
 * Returns the updated list of merged seeds, and updates nb accordingly.
 */ 
seed_t* mergeOverlappingPosSeeds(int* nb, seed_t* seeds, unsigned minOverlap);

/**
 * Merges the seeds contained in the structure of size *nb, pointed
 * by seeds if their sequences overlap on a greater length than minOverap.
 * Returns the updated list of merged seeds, and updates nb accordingly.
 */ 
seed_t* mergeOverlappingSeqSeeds(int* nb, seed_t* seeds, int minOverlap);

/**
 * Counts the number of alignments reported in the file file.
 */ 
int countAlignments(char* file);

/**
 * Returns the length of the template of identifier tpl.
 * len: length of the identifier
 */ 
int getTemplateLength(char* tpl, int len);

/**
 * Reads the alignments reported in the file file, and stores them
 * in the structure of size n_seeds, pointed by seeds.
 */
void readAlignmentFile(char* file, int n_seeds, seed_t** seeds);

/**
 * Frees the structure seeds of size size.
 */ 
void freeAlignments(int size, seed_t* seeds);

/**
 * Returns the number and list of seeds contained in the alignment
 * file file, after the two steps of processing, allowing a minimum
 * overlap of length minOverap.
 */ 
std::pair<int, seed_t*> processSeeds(char* tpl, unsigned minOverlap);


