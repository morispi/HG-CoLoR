#include "seedsLinking.h"
#include "seedsMerging.h"

namespace SLRgen {
    PgSAIndexStandard* pgsaIndex;

    typedef DefaultPgSAIndex<uint_reads_cnt_std, unsigned int
        , uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std,
            DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, 4>::Type> PgSAIndexStandardImpl;
    
    /**
     * Extends the SLR on the left.
     */ 
    void extendLeft(PgSAIndexStandard* index, int extLen, string &LR, int KLEN, int minoverlap) {
		vector<StandardOccurrence> q3res;
		vector<StandardOccurrence>::iterator it;
		int klen = KLEN;
		string kmer;
		
		// Get the leftmost k-mer and search for a path in the graph
		do {
			klen--;
			kmer = LR.substr(0, klen);
			pgsaIndex->reportOccurrences(kmer, q3res);
		} while(q3res.size() == 0 && klen > minoverlap);
		
		it = q3res.begin();
		int dist = 0;
		int read;
		int pos;
		
		// Keep traversing the graph while the template's border or a branching path aren't reached
		while (q3res.size() == 2 && it != q3res.end() && dist < extLen) {
			read = it->first;
			pos = it->second;
			if (pos == KLEN - klen) {
				string curRead = pgsaIndex->getReadVirtual(read);
				LR = curRead.substr(0, KLEN - klen) + LR;
				dist = dist + KLEN - klen;
				klen = KLEN;
				// Update current k-mer, and search for a path to continue the extension with
				do {
					klen--;
					kmer = LR.substr(0, klen);
					pgsaIndex->reportOccurrences(kmer, q3res);
				} while(q3res.size() == 0 && klen > minoverlap);
				it = q3res.begin();
			} else {
				it++;
			}
		}	
	}
	
	/**
     * Extends the SLR on the right.
     */
    void extendRight(PgSAIndexStandard* index, int extLen, string &LR, int KLEN, int minoverlap) {
		vector<StandardOccurrence> q3res;
		vector<StandardOccurrence>::iterator it;
		int klen = KLEN;
		string kmer;
		
		// Get the rightmost k-mer and search for a path in the graph
		do {
			klen--;
			kmer = LR.substr(LR.length() - klen);
			pgsaIndex->reportOccurrences(kmer, q3res);
		} while(q3res.size() == 0 && klen > minoverlap);
		
		it = q3res.begin();
		int dist = 0;
		int read;
		int pos;
		
		// Keep traversing the graph while the template's border or a branching path aren't reached
		while (q3res.size() == 2 && it != q3res.end() && dist < extLen) {
			read = it->first;
			pos = it->second;
			if (pos == 0) {
				string curRead = pgsaIndex->getReadVirtual(read);
				LR = LR + curRead.substr(klen);
				dist = dist + KLEN - klen;
				klen = KLEN;
				// Update current k-mer, and search for a path to continue the extension with
				do {
					klen--;
					kmer = LR.substr(LR.length() - klen);
					pgsaIndex->reportOccurrences(kmer, q3res);
				} while(q3res.size() == 0 && klen > minoverlap);
				it = q3res.begin();
			} else {
				it++;
			}
		}
	}
    
    /**
     * Links src to dst by traversing the graph.
     */ 
    int link(PgSAIndexStandard* index, string src, string dst, int tplLen, int klen, std::set<std::string> visited, int maxback, int* curback, int dist, string curExt, string &fRes, int KLEN, int minoverlap) {
		if (klen < minoverlap || *curback > maxback || dist > tplLen) {
				visited.clear();
				return 0;
		}
		
		vector<StandardOccurrence> q3res;
		vector<StandardOccurrence>::iterator it;
		string kmer = src.substr(src.length() - klen);
		string resPart1 = string(curExt);
		
		pgsaIndex->reportOccurrences(kmer, q3res);
		it = q3res.begin();
		int found = src == dst;	
			
		// While the destination or a braching path aren't reached, keep on 
		// traversing the graph to extend the rightmost k-mer of the source
		while (!found && q3res.size() == 2 && it != q3res.end() && dist <= tplLen) {
			int read = it->first;
			int pos = it->second;
			if (pos == 0) {
				string curRead = pgsaIndex->getReadVirtual(read);
				set<std::string>::iterator itf = visited.find(pgsaIndex->getReadVirtual(read));
				found = curRead == dst;
				if (!found && (itf == visited.end())) {
					visited.insert(pgsaIndex->getReadVirtual(read));
					kmer = curRead.substr(curRead.length() - klen);
					resPart1 = resPart1 + curRead.substr(klen, KLEN - klen);
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
		
		// If a branching path is reached, use backtracking to explore every possible path
		string tmpres;
		while (it != q3res.end() && !found && dist <= tplLen) {
			int read = it->first;
			int pos = it->second;
			if (pos == 0) {
				string curRead = pgsaIndex->getReadVirtual(read);
				set<std::string>::iterator itf = visited.find(pgsaIndex->getReadVirtual(read));
				found = curRead == dst;
				if (!found && (itf == visited.end())) {
					visited.insert(pgsaIndex->getReadVirtual(read));
					kmer = curRead.substr(curRead.length() - klen);
					tmpres = resPart1 + curRead.substr(klen, KLEN - klen);
					found = link(index, curRead, dst, tplLen, KLEN-1, visited, maxback, curback, dist + KLEN - klen, tmpres, fRes, KLEN, minoverlap);
					if (!found) {
						(*curback)++;
						visited.erase(pgsaIndex->getReadVirtual(read));
						++it;
					} else {
						visited.clear();
						return 1;
					}
				} else if (found) {
					fRes = resPart1 + dst.substr(klen);
					visited.clear();
					return 1;
				} else {
					++it;
				}
			} else {
				++it;
			}
		}
		
		// If the source couldn't be linked to the destination, try again, allowing
		// a lower overlap length
		if (!found) {
			if (klen - 1 >= minoverlap && dist < tplLen) {
				(*curback)++;
				return link(index, src, dst, tplLen, klen-1, visited, maxback, curback, dist, curExt, fRes, KLEN, minoverlap);
			} else {
				fRes = string();
				visited.clear();
				return 0;
			}
		}
		
		// This part is reached if the source could be linked to the destination without
		// encountering any branching path
		fRes = resPart1 + dst.substr(klen);
		visited.clear();
		return 1;
	}

	/**
	 * Links together the seeds contained in the file tolink, and outputs the corrected long read.
	 */ 
    void generateSLR(PgSAIndexStandard* index, int KLEN, char* tolink, string tplName, int seedsoverlap, int minoverlap, int backtracks, int seedskips) {
		pgsaIndex = index;
		int klen = KLEN-1;
		char line[4096];
		string src, dst, k1, k2;
		int fragments = 1;
		int skippedSeeds = 0;
		int firstSkippedSeed = -1;
		seed_t curSeed;
		int idSeed, posBeg, posSrc, posDst, tplLen, dist, curback, linked;
		string fRes = string();
		
		std::vector<seed_t> seeds = processSeeds(tolink, seedsoverlap);
		idSeed = 0;
		curSeed = seeds[idSeed];
		posSrc = curSeed.pos;
		posBeg = posSrc;
		tplLen = curSeed.tlen;
		src = curSeed.seq;
		std::set<std::string> visited;
		
		idSeed++;
		while (idSeed < seeds.size()) {
			curSeed = seeds[idSeed];
			posDst = curSeed.pos;
			dist = posDst - posSrc - src.length();
			dst = curSeed.seq;
			
			k1 = src.substr(src.length() - KLEN);
			k2 = dst.substr(0, KLEN);
			
			string tmpRes = string();
			curback = 0;
			linked = link(index, k1, k2, tplLen, KLEN-1, visited, backtracks, &curback, 0, src, tmpRes, KLEN, minoverlap);
			
			if (linked != 0) {
				tmpRes = tmpRes + dst.substr(KLEN);
				if (fRes.empty()) {
					fRes = tmpRes;
				} else {
					fRes = fRes + tmpRes.substr(src.length());
				}
				src = dst;
				posSrc = posDst;
			} else {
				// No seeds linked so far, skip the source
				if (fRes.empty()) {
					src = dst;
					posSrc = posDst;
					posBeg = posSrc;
				// Seeds have been linked previously, skip the source if the allowed number of skips isn't reached
				} else if (skippedSeeds < seedskips) {
					skippedSeeds++;
					if (firstSkippedSeed == -1) {
						firstSkippedSeed = idSeed;
					}
				} else {
					// Couldn't link src to dst after skipping the allowed number of seeds, so, fragment the corrected long read
					if (posBeg > 0) {
						extendLeft(index, posBeg, fRes, KLEN, minoverlap);
					}
					if (tplLen - posSrc - src.length() > 0) {
						extendRight(index, tplLen - posSrc - src.length(), fRes, KLEN, minoverlap);
					}
					cout << ">" << tplName << "_" << fragments << endl << fRes << endl;
					fragments++;
					idSeed = firstSkippedSeed;
					curSeed = seeds[idSeed];
					src = curSeed.seq;
					posSrc = curSeed.pos;
					posBeg = posSrc;
					fRes = string();
					firstSkippedSeed = -1;
					skippedSeeds = 0;
				}
			}	
			idSeed++;
		}
		
		// Multiple seeds were mapped on the template and were linked
		if (!fRes.empty()) {
			if (posBeg > 0) {
				extendLeft(index, posBeg, fRes, KLEN, minoverlap);
			}
			if (tplLen - posSrc - src.length() > 0) {
				extendRight(index, tplLen - posSrc - src.length(), fRes, KLEN, minoverlap);
			}
			if (fragments == 1) {
				cout << ">" << tplName << endl << fRes << endl;
			} else {
				cout << ">" << tplName << "_" << fragments << endl << fRes << endl;
			}
		}
		
		// Only one seed was mapped on the template
		if (seeds.size() < 2) {
			int srcLen = src.length();
			if (posBeg > 0) {
				extendLeft(index, posBeg, src, KLEN, minoverlap);
			}
			if (tplLen - posSrc - srcLen > 0) {
				extendRight(index, tplLen - posSrc - srcLen, src, KLEN, minoverlap);
			}
			cout << ">" << tplName << endl << src << endl;
		}
     }
}
