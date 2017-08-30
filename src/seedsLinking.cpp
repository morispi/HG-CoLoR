#include "seedsLinking.h"
#include "seedsMerging.h"

namespace CLRgen {
    PgSAIndexStandard* pgsaIndex;
    string tplId;
    int tplLen;
    int kMerSize;
    int minOverlap;
    int maxBackTracks;

    typedef DefaultPgSAIndex<uint_reads_cnt_std, unsigned int
        , uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std,
            DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, 4>::Type> PgSAIndexStandardImpl;
    
    
    vector<StandardOccurrence> getNeighbours(string kMer, int left) {
		vector<StandardOccurrence> neighbours;
		int read, pos, startPos;
		vector<StandardOccurrence> q3res;
		pgsaIndex->reportOccurrences(kMer, q3res);
		startPos = left == 1 ? kMerSize - kMer.size() : 0;
		vector<StandardOccurrence>::iterator it;
		it = q3res.begin();
		
		while (it != q3res.end()) {
			read = it->first;
			pos = it->second;
			if (pos == startPos) {
				neighbours.push_back(make_pair(read, pos));
			}
			it++;
		}		
		return neighbours;
	}
    
    void extendLeft(int extLen, string &LR) {
		vector<StandardOccurrence> q3res;
		vector<StandardOccurrence>::iterator it;
		int curK = kMerSize;
		int dist = 0;
		int read;
		int pos;
		string curRead;
		string kMer;
		
		// Get the leftmost k-mer and search for a path in the graph
		do {
			curK--;
			kMer = LR.substr(0, curK);
			q3res = getNeighbours(kMer, 1);
		} while(q3res.size() == 0 && curK > minOverlap);
		
		it = q3res.begin();
		
		
		// Keep traversing the graph while the template's border or a branching path aren't reached
		while (q3res.size() == 1 && it != q3res.end() && dist < extLen) {
			read = it->first;
			curRead = pgsaIndex->getReadVirtual(read);
			LR = curRead.substr(0, kMerSize - curK) + LR;
			dist = dist + kMerSize - curK;
			curK = kMerSize;
			// Update current k-mer, and search for a path to continue the extension with
			do {
				curK--;
				kMer = LR.substr(0, curK);
				q3res = getNeighbours(kMer, 1);
			} while(q3res.size() == 0 && curK > minOverlap);
			it = q3res.begin();
		}	
	}
	
	void extendRight(int extLen, string &LR) {
		vector<StandardOccurrence> q3res;
		vector<StandardOccurrence>::iterator it;
		int curK = kMerSize;
		int dist = 0;
		int read;
		int pos;
		string curRead;
		string kMer;
		
		// Get the rightmost k-mer and search for a path in the graph
		do {
			curK--;
			kMer = LR.substr(LR.length() - curK);
			q3res = getNeighbours(kMer, 0);
		} while(q3res.size() == 0 && curK > minOverlap);
		
		it = q3res.begin();
		
		// Keep traversing the graph while the template's border or a branching path aren't reached
		while (q3res.size() == 1 && it != q3res.end() && dist < extLen) {
			read = it->first;
			curRead = pgsaIndex->getReadVirtual(read);
			LR = LR + curRead.substr(curK);
			dist = dist + kMerSize - curK;
			// Update current k-mer, and search for a path to continue the extension with
			curK = kMerSize;
			do {
				curK--;
				kMer = LR.substr(LR.length() - curK);
				q3res = getNeighbours(kMer, 0);
			} while(q3res.size() == 0 && curK > minOverlap);
			it = q3res.begin();
		}
	}
	
	int link(string src, string dst, int curK, set<int> &visited, int* curback, int dist, string curExt, string &fullLinkingSeq) {
		if (curK < minOverlap || *curback > maxBackTracks || dist > tplLen) {
				return 0;
		}
		
		vector<StandardOccurrence> q3res;
		vector<StandardOccurrence>::iterator it;
		int read;
		int found = src == dst;
		string curRead;
		string kMer;
		string resPart1 = string(curExt);
		string tmpres;
		set<int>::iterator itf;
		
		// Get the rightmost k-mer and search for a path in the graph
		do {
			curK--;
			kMer = src.substr(src.length() - curK);
			q3res = getNeighbours(kMer, 0);
		} while(!found && q3res.size() == 0 && curK > minOverlap);
		
		it = q3res.begin();
			
		/*
		 * While the destination or a braching path aren't reached, keep on 
		 * traversing the graph to extend the rightmost k-mer of the source
		 */
		while (!found && q3res.size() == 1 && it != q3res.end() && dist <= tplLen) {
			read = it->first;
			curRead = pgsaIndex->getReadVirtual(read);
			itf = visited.find(read);
			found = curRead == dst;
			if (!found && (itf == visited.end())) {
				visited.insert(read);
				resPart1 = resPart1 + curRead.substr(curK, kMerSize - curK);
				dist = dist + kMerSize - curK;

				// Update current k-mer, and search for a path to continue the extension with
				curK = kMerSize;
				do {
					curK--;
					kMer = curRead.substr(curRead.length() - curK);
					q3res = getNeighbours(kMer, 0);
				} while(q3res.size() == 0 && curK > minOverlap);
				
				it = q3res.begin();
			} else {
				it++;
			}
		}
		
		/*
		 * If a branching path is reached, explore the different possible
		 * paths with backtracking
		 */
		while (!found && q3res.size() > 1 && it != q3res.end() && dist <= tplLen) {
			read = it->first;
			curRead = pgsaIndex->getReadVirtual(read);
			itf = visited.find(read);
			found = curRead == dst;
			if (!found && (itf == visited.end())) {
				visited.insert(read);
				kMer = curRead.substr(curRead.length() - curK);
				tmpres = resPart1 + curRead.substr(curK, kMerSize - curK);
				found = link(curRead, dst, kMerSize, visited, curback, dist + kMerSize - curK, tmpres, fullLinkingSeq);
				if (!found) {
					(*curback)++;
					++it;
				} else {
					return 1;
				}
			} else if (found) {
				fullLinkingSeq = resPart1 + dst.substr(curK);
				return 1;
			} else {
				++it;
			}
		}
		
		/*
		 * If the source couldn't be linked to the destination, try again with a
		 * smaller overlap length
		 */
		if (!found) {
			if (curK - 1 >= minOverlap && dist < tplLen) {
				(*curback)++;
				return link(src, dst, curK-1, visited, curback, dist, curExt, fullLinkingSeq);
			} else {
				fullLinkingSeq = string();
				return 0;
			}
		}
		
		/*
		 * This part is reached if the source could be linked to the destination without
		 * encountering any branching path
		 */
		fullLinkingSeq = resPart1 + dst.substr(curK);
		return 1;
	}

	/**
	 * Links together the seeds contained in the file tolink, and outputs the corrected long read.
	 */ 
    void generateCLR(PgSAIndexStandard* index, int k, string tolink, string templateId, int seedsoverlap, int minoverlap, int backtracks, int seedskips) {
		// global variables
		pgsaIndex = index;
		tplId = templateId;
		minOverlap = minoverlap;
		kMerSize = k;
		maxBackTracks = backtracks;
		
		int curK = k - 1;
		int fragments = 1;
		int skippedSeeds = 0;
		int firstSkippedSeed = -1;
		int idSeed, posBeg, posSrc, posDst, dist, curback, linked;
		seed_t curSeed;
		string src, dst, kMer1, kMer2, fullLinkingSeq;
		ostringstream fRes;
		set<int> visited;
	
		
		vector<seed_t> seeds = processSeeds(tolink, seedsoverlap);
		idSeed = 0;
		curSeed = seeds[idSeed];
		posSrc = curSeed.pos;
		posBeg = posSrc;
		tplLen = curSeed.tlen;
		src = curSeed.seq;
		string curLinkingSeq;
		
		if (seedskips > (seeds.size() - 1) - idSeed - 1) {
			seedskips = (seeds.size() - 1) - idSeed - 1;
		}
		
		idSeed++;
		// Iterate through the seeds and link them
		while (idSeed < seeds.size()) {
			curSeed = seeds[idSeed];
			posDst = curSeed.pos;
			dist = posDst - posSrc - src.length();
			dst = curSeed.seq;
			
			kMer1 = src.substr(src.length() - kMerSize);
			kMer2 = dst.substr(0, kMerSize);
			
			curLinkingSeq = string();
			curback = 0;
			linked = link(kMer1, kMer2, kMerSize, visited, &curback, 0, src, curLinkingSeq);
			visited.clear();
			// Seeds could be linked, update the result string
			if (linked != 0) {
				curLinkingSeq = curLinkingSeq + dst.substr(kMerSize);
				if (fullLinkingSeq.empty()) {
					fullLinkingSeq = curLinkingSeq;
				} else {
					fullLinkingSeq = fullLinkingSeq + curLinkingSeq.substr(src.length());
				}
				src = dst;
				posSrc = posDst;
				firstSkippedSeed = -1;
				skippedSeeds = 0;
				if (seedskips > (seeds.size() - 1) - idSeed - 1) {
					seedskips = (seeds.size() - 1) - idSeed - 1;
				}
			// Seeds couldn't be linked, skip a seed or fragment the correct long read 
			} else {
				// No seeds linked so far, skip the source
				if (fullLinkingSeq.empty()) {
					src = dst;
					posSrc = posDst;
					posBeg = posSrc;
					if (seedskips > (seeds.size() - 1) - idSeed - 1) {
						seedskips = (seeds.size() - 1) - idSeed - 1;
					}
				// Seeds have been linked previouslyskip the destination if the allowed number of skips isn't reached
				} else if (skippedSeeds < seedskips) {
					skippedSeeds++;
					if (firstSkippedSeed == -1) {
						firstSkippedSeed = idSeed;
					}
				} else {
					// Couldn't link src to dst after skipping the allowed number of seeds, so fragment the corrected long read
					if (posBeg > 0) {
						extendLeft(posBeg, fullLinkingSeq);
					}
					if (tplLen - posSrc - src.length() > 0) {
						extendRight(tplLen - posSrc - src.length(), fullLinkingSeq);
					}
					fRes << ">" << tplId << "_" << fragments << endl << fullLinkingSeq << endl;
					fragments++;
					if (firstSkippedSeed != -1) {
						idSeed = firstSkippedSeed;
					}
					curSeed = seeds[idSeed];
					src = curSeed.seq;
					posSrc = curSeed.pos;
					posBeg = posSrc;
					fullLinkingSeq = string();
					firstSkippedSeed = -1;
					skippedSeeds = 0;
					if (seedskips > (seeds.size() - 1) - idSeed - 1) {
						seedskips = (seeds.size() - 1) - idSeed - 1;
					}
				}
			}
			idSeed++;
		}

		// Multiple seeds were mapped on the template and were linked
		if (!fullLinkingSeq.empty()) {
			if (posBeg > 0) {
				extendLeft(posBeg, fullLinkingSeq);
			}
			if (tplLen - posSrc - src.length() > 0) {
				extendRight(tplLen - posSrc - src.length(), fullLinkingSeq);
			}
			if (fragments == 1) {
				fRes << ">" << tplId << endl << fullLinkingSeq << endl;
			} else {
				fRes << ">" << tplId << "_" << fragments << endl << fullLinkingSeq << endl;
			}
		}

		// Only one seed was mapped on the template, extend it
		if (seeds.size() < 2) {
			int srcLen = src.length();
			if (posBeg > 0) {
				extendLeft(posBeg, src);
			}
			if (tplLen - posSrc - srcLen > 0) {
				extendRight(tplLen - posSrc - srcLen, src);
			}
			//~ cout << ">" << tplId << endl << src << endl;
			fRes << ">" << tplId << endl << src << endl;
			//~ fRes = fRes + ">" + tplId + endl + src + endl;
		}
		
		cout << fRes.str();
     }
     
     
     
}
