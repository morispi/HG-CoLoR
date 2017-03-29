#include "seedsLinking.h"
#include "seedsMerging.h"

namespace PgSATest {
    
    int repeat = 11;
    int testKmersCount = 1;

    int kValue;
    
    PgSAIndexStandard* pgsaIndex;
    string* testkmers = 0;
    std::pair<t_reads_c, t_read>* testFactors;
	
	bool sortByOccurrencesNb(const pair<StandardOccurrence, unsigned int>& p1, const pair<StandardOccurrence, unsigned int>& p2) {
		return p1.second > p2.second;
	}
    

    typedef DefaultPgSAIndex<uint_reads_cnt_std, unsigned int
        , uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std,
            DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, 4>::Type> PgSAIndexStandardImpl;

    void doFilterTTTs() {
        int tttsCount = 0;
        for (int i = 0; i < testKmersCount; i++) {
            if (testkmers[i].find_first_of("ACGN") == string::npos) {
                tttsCount++;
                int idx = (i + testKmersCount - 1) % testKmersCount;
                testkmers[i] = testkmers[idx];
                testFactors[i] = testFactors[idx];
            }
        }
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void drawDataK(ReadsSetInterface<uint_read_len, uint_reads_cnt>* readsSet, unsigned short k) {
        if (testkmers != 0)
            delete[]testkmers;
        if (testFactors != 0)
            delete[]testFactors;
            
        testFactors = generateTagFactors(readsSet->readsCountVirtual(), readsSet->maxReadLengthVirtual(), k, testKmersCount);
        testkmers = generateTextFactors(readsSet, k, testFactors, testKmersCount);
        
        kValue = k;

    }

    inline const string getKmer(t_reads_c i) {
        const std::pair<t_reads_c, t_read> pair = testFactors[i];
        return pgsaIndex->getReadVirtual(pair.first).substr(pair.second, kValue);
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    void drawDataK(ReadsSetInterface<uint_read_len, uint_reads_cnt>* readsSet, unsigned short k, bool pairScrambled ) {
        if (testkmers != 0)
            delete[]testkmers;
        if (testFactors != 0)
            delete[]testFactors;
            
        if (pairScrambled) 
            testFactors = generatePairScrambledTagFactors(readsSet->readsCountVirtual(), readsSet->maxReadLengthVirtual(), k, testKmersCount);
        else        
            testFactors = generateTagFactors(readsSet->readsCountVirtual(), readsSet->maxReadLengthVirtual(), k, testKmersCount); 
            
        testkmers = generateTextFactors(readsSet, k, testFactors, testKmersCount);
        
        kValue = k;   
    }
    
    void extendLeft(PgSAIndexStandard* index, bool pairScrambled, int extLen, char** LR, int KLEN, int minoverlap) {
		vector<StandardOccurrence> q3res;
		vector<StandardOccurrence>::iterator it;
		int klen = KLEN;
		char kmer[KLEN];
		
		// Search for overlapping k-mers
		do {
			klen--;
			strncpy(kmer, *LR, klen);
			kmer[klen] = '\0';
			drawDataK(index, klen, pairScrambled);
			pgsaIndex->reportOccurrences(kmer, q3res);
		} while(q3res.size() == 0 && klen > minoverlap);
		
		it = q3res.begin();
		int dist = 0;
		int read;
		int pos;
		
		// Extend while the template's border or an ambiguity aren't reached
		while (q3res.size() == 2 && it != q3res.end() && dist < extLen) {
			read = it->first;
			pos = it->second;
			if (pos == KLEN - klen) {
				char curRead[pgsaIndex->getReadVirtual(read).length() + 1];
				pgsaIndex->getReadVirtual(read).copy(curRead, pgsaIndex->getReadVirtual(read).length(), 0);
				curRead[pgsaIndex->getReadVirtual(read).length()] = '\0';
				char tmp[strlen(*LR) + 1];
				strncpy(tmp, *LR, strlen(*LR) + 1);
				free(*LR);
				*LR = (char*) malloc(strlen(tmp) + KLEN - klen + 1);
				strncpy(*LR, curRead, KLEN - klen);
				(*LR)[KLEN - klen] = '\0';
				strncat(*LR, tmp, strlen(tmp));
				dist = dist + KLEN - klen;
				klen = KLEN;
				// Update current k-mer, and search for the next overlapping k-mers
				do {
					klen--;
					strncpy(kmer, *LR, klen);
					kmer[klen] = '\0';
					drawDataK(index, klen, pairScrambled);
					pgsaIndex->reportOccurrences(kmer, q3res);
				} while(q3res.size() == 0 && klen > minoverlap);
				it = q3res.begin();
			} else {
				it++;
			}
		}	
	}
	
    void extendRight(PgSAIndexStandard* index, bool pairScrambled, int extLen, char** str, int KLEN, int minoverlap) {
		vector<StandardOccurrence> q3res;
		vector<StandardOccurrence>::iterator it;
		int klen = KLEN;
		char kmer[KLEN];
		
		// Search for overlapping k-mers
		do {
			klen--;
			strncpy(kmer, *str + strlen(*str) - klen, klen + 1);
			drawDataK(index, klen, pairScrambled);
			pgsaIndex->reportOccurrences(kmer, q3res);
		} while(q3res.size() == 0 && klen > minoverlap);
		
		it = q3res.begin();
		int dist = 0;
		int read;
		int pos;
		
		// Extend while the template's border or an ambiguity aren't reached
		while (q3res.size() == 2 && it != q3res.end() && dist < extLen) {
			read = it->first;
			pos = it->second;
			if (pos == 0) {
				char curRead[pgsaIndex->getReadVirtual(read).length() + 1];
				pgsaIndex->getReadVirtual(read).copy(curRead, pgsaIndex->getReadVirtual(read).length(), 0);
				curRead[pgsaIndex->getReadVirtual(read).length()] = '\0';
				*str = (char*) realloc(*str, strlen(*str) + KLEN - klen + 1);
				strncat(*str, curRead + klen, KLEN - klen);
				dist = dist + KLEN - klen;
				klen = KLEN;
				// Update current k-mer, and search for the next overlapping k-mers
				do {
					klen--;
					strncpy(kmer, *str + strlen(*str) - klen, klen + 1);
					drawDataK(index, klen, pairScrambled);
					pgsaIndex->reportOccurrences(kmer, q3res);
				} while(q3res.size() == 0 && klen > minoverlap);
				it = q3res.begin();
			} else {
				it++;
			}
		}
	}
    
    int link(PgSAIndexStandard* index, char* src, char* dst, int gapSize, int klen, bool pairScrambled, std::set<std::string> visited, int maxback, int* curback, int dist, char* res, char** fres, int KLEN, int minoverlap) {
		if (klen < minoverlap || *curback > maxback || dist > gapSize) {
				visited.clear();
				return 0;
		}
		
		drawDataK(index, klen, pairScrambled);
		vector<StandardOccurrence> q3res;
		vector<StandardOccurrence>::iterator it;
		char kmer[klen+1];
		char* resPart1 = (char*) malloc(strlen(res) + 1);
		strncpy(resPart1, res, strlen(res) + 1);
		strncpy(kmer, src + strlen(src) - klen, klen + 1);
		pgsaIndex->reportOccurrences(kmer, q3res);
		it = q3res.begin();
		int found = strcmp(src, dst) == 0 ? 1 : 0;
		
		// While there's no ambiguity, we keep on extending the rightmost k-mer of the left seed
		while (!found && q3res.size() == 2 && it != q3res.end() && dist <= gapSize) {
			int read = it->first;
			int pos = it->second;
			if (pos == 0) {
				char curRead[pgsaIndex->getReadVirtual(read).length() + 1];
				pgsaIndex->getReadVirtual(read).copy(curRead, pgsaIndex->getReadVirtual(read).length(), 0);
				curRead[pgsaIndex->getReadVirtual(read).length()] = '\0';
				set<std::string>::iterator itf = visited.find(pgsaIndex->getReadVirtual(read));
				found = strcmp(dst, curRead) == 0 ? 1 : 0;
				if (!found && (itf == visited.end())) {
					visited.insert(pgsaIndex->getReadVirtual(read));
					char src2[strlen(curRead) + 1];
					strncpy(src2, curRead, strlen(curRead) + 1);
					strncpy(kmer, src2 + strlen(src2) - klen, klen + 1);
					resPart1 = (char*) realloc(resPart1, strlen(resPart1) + KLEN - klen + 1);
					strncat(resPart1, curRead + klen, KLEN - klen);
					pgsaIndex->reportOccurrences(kmer, q3res);
					it = q3res.begin();
					dist = dist + KLEN - klen;
				} else {
					it++;
				}
			} else {
				++it;
			}
		}
		
		// If multiple k-mers allow the extension of the current k-mer, we use backtracking
		char* tmpres;
		while (it != q3res.end() && !found && dist <= gapSize) {
			int read = it->first;
			int pos = it->second;
			if (pos == 0) {
				char curRead[pgsaIndex->getReadVirtual(read).length() + 1];
				pgsaIndex->getReadVirtual(read).copy(curRead, pgsaIndex->getReadVirtual(read).length(), 0);
				curRead[pgsaIndex->getReadVirtual(read).length()] = '\0';
				set<std::string>::iterator itf = visited.find(pgsaIndex->getReadVirtual(read));
				found = strcmp(dst, curRead) == 0 ? 1 : 0;
				if (!found && (itf == visited.end())) {
					visited.insert(pgsaIndex->getReadVirtual(read));
					char src2[strlen(curRead) + 1];
					strncpy(src2, curRead, strlen(curRead) + 1);
					strncpy(kmer, src2 + strlen(src2) - klen, klen + 1);
					tmpres = (char*) malloc(strlen(resPart1) + KLEN - klen + 1);
					strncpy(tmpres, resPart1, strlen(resPart1) + 1);
					strncat(tmpres, curRead + klen, KLEN - klen);
					found = link(index, src2, dst, gapSize, KLEN-1, pairScrambled, visited, maxback, curback, dist + KLEN - klen, tmpres, fres, KLEN, minoverlap);
					free(tmpres);
					if (!found) {
						(*curback)++;
						visited.erase(pgsaIndex->getReadVirtual(read));
						++it;
					} else {
						visited.clear();
						return 1;
					}
				} else if (found) {
				    *fres = (char*) malloc(strlen(resPart1) + strlen(dst) - klen + 1);
					strncpy(*fres, resPart1, strlen(resPart1) + 1);
					strncat(*fres, dst + klen, strlen(dst) - klen);
					visited.clear();
					return 1;
				} else {
					++it;
				}
			} else {
				++it;
			}
		}
		
		// If seeds couldn't be linked, we try again with a lower overlap length between k-mers
		if (!found) {
			if (klen - 1 >= minoverlap && dist < gapSize) {
				(*curback)++;
				return link(index, src, dst, gapSize, klen-1, pairScrambled, visited, maxback, curback, dist, res, fres, KLEN, minoverlap);
			} else {
				*fres = NULL;
				visited.clear();
				return 0;
			}
		}
		
		// We reach this part if the two seeds could be linked with no ambiguity
		*fres = (char*) malloc(strlen(resPart1) + strlen(dst) - klen + 1);
		strncpy(*fres, resPart1, strlen(resPart1) + 1);
		strncat(*fres, dst + klen, strlen(dst) - klen);
		visited.clear();
		return 1;
	}

    void runTest(PgSAIndexStandard* index, int _repeat, int _testKmersCount, unsigned short k, bool pairScrambled, bool filterTTTs, bool byPosition, int KLEN, char* tolink, char* tplName, int seedsoverlap, int minoverlap, int backtracks, int seedskips) {
		pgsaIndex = index;
        testKmersCount = _testKmersCount;
		int klen = KLEN-1;
        drawDataK(index, klen, pairScrambled);
	
		int skippedSeeds = 0;
		char line[4096];
		char* src = NULL;
		char k1[KLEN+1];
		char k2[KLEN+1];
		int gapSize;
		int dist;
		int posBeg;
		int fragments = 1;
		int firstSkippedSeed = -1;
		std::pair<int, seed_t*> processedSeeds = processSeeds(tolink, seedsoverlap);
		int nbSeeds = processedSeeds.first;
		seed_t* seeds = processedSeeds.second;
		
		int idSeed = 0;
		seed_t curSeed = *seeds;
		int posTpl = curSeed.pos;
		posBeg = posTpl;
		gapSize = curSeed.tlen;
		src = (char*) malloc(strlen(curSeed.seq) + 1);
		strncpy(src, strToUpper(curSeed.seq, strlen(curSeed.seq)), strlen(curSeed.seq) + 1);
		
		
		int outSrc = strlen(src);
		std::set<std::string> visited;
		char* fRes = NULL;
		
		idSeed++;
		while (idSeed < nbSeeds) {
			outSrc = 0;
			curSeed = *(seeds + idSeed);
			int tmpposTpl = curSeed.pos;
			dist = tmpposTpl - posTpl - strlen(src);
			char dst[strlen(curSeed.seq) + 1];
			strncpy(dst, strToUpper(curSeed.seq, strlen(curSeed.seq)), strlen(curSeed.seq) + 1);
			
			
			strncpy(k1, src + strlen(src) - KLEN, KLEN + 1);
			strncpy(k2, dst, KLEN);
			k2[KLEN] = '\0';
			
			char* src2 = (char*) malloc(strlen(src) + 1);;
			strncpy(src2, src, strlen(src) + 1);
			

			char* res2 = NULL;
			int klen;
			int curback = 0;
			int linked = link(index, k1, k2, gapSize, KLEN-1, pairScrambled, visited, backtracks, &curback, 0, src2, &res2, KLEN, minoverlap);
			if (linked != 0) {
				res2 = (char*) realloc(res2, strlen(res2) + strlen(dst) - KLEN + 1);
				strncat(res2, dst + KLEN, strlen(dst) - KLEN);
				if (fRes == NULL) {
					fRes = (char*) malloc(strlen(res2) + 1);
					strncpy(fRes, res2, strlen(res2) + 1);
				} else {
					fRes = (char*) realloc(fRes, strlen(fRes) + strlen(res2) - strlen(src) + 1);
					strncat(fRes, res2 + strlen(src), strlen(res2) - strlen(src));
				}
				free(src);
				src = (char*) malloc(strlen(dst) + 1);
				strncpy(src, dst, strlen(dst) + 1);
				posTpl = tmpposTpl;
			} else {
				if (fRes == NULL) {
					free(src);
					src = (char*) malloc(strlen(dst) + 1);
					strncpy(src, dst, strlen(dst) + 1);
					posTpl = tmpposTpl;
					posBeg = posTpl;
				} else if (skippedSeeds < seedskips) {
					skippedSeeds++;
					if (firstSkippedSeed == -1) {
						firstSkippedSeed = idSeed;
					}
				} else {
					// extend, output, fragment
					if (posBeg > 0) {
						extendLeft(index, pairScrambled, posBeg, &fRes, KLEN, minoverlap);
					}
					if (gapSize - posTpl - strlen(src) > 0) {
						extendRight(index, pairScrambled, gapSize - posTpl - strlen(src), &fRes, KLEN, minoverlap);
					}
					cerr << ">" << tplName << "_" << fragments << endl;
					cerr << fRes << endl;
					fragments++;
					idSeed = firstSkippedSeed;
					curSeed = *(seeds + idSeed);
					free(src);
					src = (char*) malloc(strlen(curSeed.seq) + 1);
					strncpy(src, strToUpper(curSeed.seq, strlen(curSeed.seq)), strlen(curSeed.seq) + 1);
					posTpl = curSeed.pos;
					posBeg = posTpl;
					free(fRes);
					fRes = NULL;
					firstSkippedSeed = -1;
					skippedSeeds = 0;
				}
			}
			
			free(src2);
			src2 = NULL;
			free(res2);
			res2 = NULL;
			klen = 0;
			idSeed++;
		}
		freeAlignments(nbSeeds, seeds);
		
		if (fRes != NULL) {
			if (posBeg > 0) {
				extendLeft(index, pairScrambled, posBeg, &fRes, KLEN, minoverlap);
			}
			if (gapSize - posTpl - strlen(src) > 0) {
				extendRight(index, pairScrambled, gapSize - posTpl - strlen(src), &fRes, KLEN, minoverlap);
			}
			if (fragments == 1) {
				cerr << ">" << tplName << endl;
			} else {
				cerr << ">" << tplName << "_" << fragments << endl;
			}
			cerr << fRes << endl;
		}
		
		if (outSrc != 0) {
			if (posBeg > 0) {
				extendLeft(index, pairScrambled, posBeg, &src, KLEN, minoverlap);
			}
			if (gapSize - posTpl - outSrc > 0) {
				extendRight(index, pairScrambled, gapSize - posTpl - outSrc, &src, KLEN, minoverlap);
			}
			cerr << ">" << tplName << endl;
			cerr << src << endl;
		}
		
		free(src);
		free(fRes);
		fRes = NULL;
		
        delete[]testkmers;
        testkmers = 0;
    }

    
}
