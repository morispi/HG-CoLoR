#include "seedsLinking.h"
#include "seedsMerging.h"
#include "kmc_query/kmc_query.h"
#include <chrono>
#include <mutex>
#include <future>

namespace CLRgen {
    PgSAIndexStandard* pgsaIndex;
    unsigned maxOrder;
    unsigned minOrder;
    unsigned seedsOverlap;
    unsigned maxBranches;
    recursive_mutex occMtx;
    mutex outMtx, PgSAMtx;
    string tmpDir;
    unsigned maxSeedsSkips;

    typedef DefaultPgSAIndex<uint_reads_cnt_std, unsigned int
        , uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std,
            DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, 4>::Type> PgSAIndexStandardImpl;
    
    // Returns the number of difference between s1 and s2
    int getDifferences(string s1, string s2) {
		int diff = 0;
		for (unsigned i = 0; i < s1.size(); i++) {
			if (s1[i] != s2[i]) {
				diff++;
			}
		}
		return diff;
	}
    
    vector<string> getNeighbours(string kMer, int left) {
		set<string> neighbours;
		int read, pos;
		vector<StandardOccurrence> qRes;
		string f = left == 0 ? kMer.substr(1) : kMer.substr(0,kMer.length() - 1);
		PgSAMtx.lock();
		pgsaIndex->reportOccurrences(f, qRes);
		PgSAMtx.unlock();
		vector<StandardOccurrence>::iterator it;
		it = qRes.begin();
		
		while (it != qRes.end()) {
			read = it->first;
			pos = it->second;
			if (left == 0 && pos + f.length() < maxOrder) {
				neighbours.insert(pgsaIndex->getReadVirtual(read).substr(pos));
			} else if (left == 1 && pos - 1 >= 0) {
				neighbours.insert(pgsaIndex->getReadVirtual(read).substr(0, pos + f.length()));
			}
			it++;
		}
		
		vector<string> fres;
		set<string>::iterator n, i;
		n = neighbours.begin();
		i = neighbours.begin();
		while (n != neighbours.end()) {
			i++;
			while (i != neighbours.end() && i->find(*n) != std::string::npos) {
				n++;
				i++;
			}
			fres.push_back(*n);
			n++;
		}
		
		return fres;
	}
    
    void extendLeft(unsigned extLen, string &LR) {
		vector<string> neighbours;
		vector<string>::iterator it;
		unsigned curK = maxOrder;
		unsigned dist = 0;

		// Get the leftmost k-mer and search for a path in the graph
		neighbours = getNeighbours(LR.substr(0, curK), 1);
		while (curK > minOrder && neighbours.size() == 0) {
				curK--;
				neighbours = getNeighbours(LR.substr(0, curK), 1);
		}
		it = neighbours.begin();

		// Keep traversing the graph while the long reads's border or a branching path aren't reached
		while (neighbours.size() == 1 && it != neighbours.end() && dist < extLen) {
			LR = (*it).substr(0, it->length() - (curK - 1)) + LR;
			dist = dist + it->length() - (curK - 1);
			// Get the leftmost k-mer and search for a path in the graph
			curK = maxOrder;
			neighbours = getNeighbours(LR.substr(0, curK), 1);
			while (curK > minOrder && neighbours.size() == 0) {
					curK--;
					neighbours = getNeighbours(LR.substr(0, curK), 1);
			}
			it = neighbours.begin();	
		}
	}

	void extendRight(unsigned extLen, string &LR) {
		vector<string> neighbours;
		vector<string>::iterator it;
		unsigned curK = maxOrder;
		unsigned dist = 0;

		// Get the leftmost k-mer and search for a path in the graph
		neighbours = getNeighbours(LR.substr(LR.length() - curK), 0);
		while (curK > minOrder && neighbours.size() == 0) {
				curK--;
				neighbours = getNeighbours(LR.substr(LR.length() - curK), 0);
		}
		it = neighbours.begin();
		
		// Keep traversing the graph while the long reads's border or a branching path aren't reached
		while (neighbours.size() == 1 && it != neighbours.end() && dist < extLen) {
			LR = LR + (*it).substr(curK - 1);
			dist = dist + it->length() - (curK - 1);
			// Get the leftmost k-mer and search for a path in the graph
			curK = maxOrder;
			neighbours = getNeighbours(LR.substr(LR.length() - curK), 0);
			while (curK > minOrder && neighbours.size() == 0) {
					curK--;
					neighbours = getNeighbours(LR.substr(LR.length() - curK), 0);
			}
			it = neighbours.begin();
		}
	}
	
	int link(string srcSeed, string tgtSeed, unsigned curK, set<string> &visited, unsigned* curBranches, unsigned dist, string curExt, string &missingPart, unsigned LRLen) {
		if (curK <= minOrder || *curBranches > maxBranches || dist > LRLen) {
				return 0;
		}
		
		string srcAnchor = curExt.substr(curExt.length() - curK);
		string tgtAnchor = tgtSeed.substr(0, curK);
		vector<string> neighbours;
		vector<string>::iterator it;
		int found = getDifferences(srcAnchor, tgtAnchor) <= 3;
		string curRead;
		string resPart1 = string(curExt);
		set<string>::iterator itf;
		
		// Search for a path in the graph starting from the source's anchor
		neighbours = getNeighbours(srcAnchor.substr(srcAnchor.length() - curK), 0);
		while (!found && curK > minOrder && neighbours.size() == 0) {
			curK--;
			srcAnchor = curExt.substr(curExt.length() - curK);
			neighbours = getNeighbours(srcAnchor.substr(srcAnchor.length() - curK), 0);
		}
		it = neighbours.begin();

		// While the destination or a braching path aren't reached, keep on traversing the graph
		while (!found && neighbours.size() == 1 && it != neighbours.end() && dist <= LRLen) {
			curRead = *it;
			itf = visited.find(curRead);
			tgtAnchor = tgtSeed.substr(0, curRead.length());
			found = getDifferences(curRead, tgtAnchor) <= 3;
			if (!found && (itf == visited.end())) {
				visited.insert(curRead);
				resPart1 = resPart1 + curRead.substr(curK - 1);
				dist = dist + curRead.length() - (curK - 1);

				// Update the current k-mer, and search for a path in the graph
				curK = maxOrder;
				srcAnchor = resPart1.substr(resPart1.length() - curK);
				neighbours = getNeighbours(srcAnchor.substr(srcAnchor.length() - curK), 0);
				while (!found && curK > minOrder && neighbours.size() == 0) {
						curK--;
						srcAnchor = resPart1.substr(resPart1.length() - curK);
						neighbours = getNeighbours(srcAnchor.substr(srcAnchor.length() - curK), 0);
				}
				it = neighbours.begin();
			} else {
				it++;
			}
		}

		// If a branching path is reached, sort the edges according to the number of occurrences of the k-mers they lead to.
		if (!found && neighbours.size() > 1) {
			sort( neighbours.begin( ), neighbours.end( ), [ ]( string& k1, string& k2) {
				PgSAMtx.lock();
				int occ1 = k1.size() < maxOrder ? pgsaIndex->countOccurrences(k1) : getOccNb(k1);
				int occ2 = k2.size() < maxOrder ? pgsaIndex->countOccurrences(k2) : getOccNb(k2);
				PgSAMtx.unlock();
				return  occ1 > occ2;
			});
			it = neighbours.begin();
		}
		
		// If a branching path is reached, explore the different possible paths with backtracking
		while (!found && neighbours.size() > 1 && it != neighbours.end() && dist <= LRLen) {
			curRead = *it;
			itf = visited.find(curRead);
			tgtAnchor = tgtSeed.substr(0, curRead.length());
			found = getDifferences(curRead, tgtAnchor) <= 3;
			if (!found && (itf == visited.end())) {
				visited.insert(curRead);
				(*curBranches)++;
				found = link(srcSeed, tgtSeed, maxOrder, visited, curBranches, dist + curRead.length() - (curK - 1), resPart1 + curRead.substr(curK - 1), missingPart, LRLen);
				if (!found) {
					++it;
				} else {
					return 1;
				}
			} else if (found) {
				missingPart = resPart1 + tgtSeed.substr(curK - 1);
				return 1;
			} else {
				++it;
			}
		}
		
		// If the source couldn't be linked to the destination, try again with a graph of smaller order, otherwhise update the missing part and return
		if (!found) {
			if (curK - 1 > minOrder && dist < LRLen) {
				return link(srcSeed, tgtSeed, curK - 1, visited, curBranches, dist, curExt, missingPart, LRLen);
			} else {
				missingPart = string();
				return 0;
			}
		} else {
			missingPart = resPart1 + tgtSeed.substr(curK - 1);
			return 1;
		}
	}

    void generateCLRs(vector<string>& longReads) {
		unsigned fragments, skippedSeeds, posBeg, posSrc, posTgt, dist,
				 curBranches, LRLen, idSeed, seedsSkips, curLink = 0;
		int firstSkippedSeed, linked;
		string LRId, src, tgt, correctedLR;
		set<string> visited;
		vector<seed_t> seeds;
		seed_t curSeed;
		
		// Iterate through the long reads
		while (curLink < longReads.size()) {
			LRId = longReads[curLink];
			fragments = 1;
			skippedSeeds = 0;
			firstSkippedSeed = -1;
			idSeed = 0;
			posBeg = 0;
			posSrc = 0;
			posTgt = 0;
			dist = 0;
			curBranches = 0;
			linked = 0;
			LRLen = 0;
			src = "";
			tgt = "";
			correctedLR = "";
			ostringstream fRes;
			seedsSkips = maxSeedsSkips;
			
			// Merge the overlapping seeds of the current long read
			seeds = processSeeds(tmpDir + "/Alignments/" + LRId, seedsOverlap);
			
			if (seeds.size() > 0) {
				
				if (seedsSkips > (seeds.size() - 1) - idSeed - 1) {
					seedsSkips = (seeds.size() - 1) - idSeed - 1;
				}
				
				idSeed = 0;
				curSeed = seeds[idSeed];
				posSrc = curSeed.pos;
				posBeg = posSrc;
				LRLen = curSeed.tlen;
				src = curSeed.seq;
				string missingPart;
				idSeed++;
				int nextReachable;
			
				// Iterate through the seeds and link them
				while (idSeed < seeds.size()) {
					nextReachable = -1;
					curSeed = seeds[idSeed];
					posTgt = curSeed.pos;
					dist = posTgt - posSrc - src.length();
					tgt = curSeed.seq;					
					missingPart = string();
					curBranches = 0;
					linked = link(src, tgt, maxOrder, visited, &curBranches, 0, src, missingPart, LRLen);
					visited.clear();
					unsigned idTmp = idSeed + 1;
					unsigned tmpSkips = 0;
					unsigned isLinkable = idTmp == seeds.size();
					string tmpSeed;
					string tmpMissingPart = string();
					// If the target seed could be reached, check if it can be linked
					while (linked && !isLinkable && idTmp < seeds.size() && tmpSkips <= maxSeedsSkips) {
						tmpSeed = seeds[idTmp].seq;
						curBranches = 0;
						isLinkable = link(tgt, tmpSeed, maxOrder, visited, &curBranches, 0, tgt, tmpMissingPart, LRLen);
						visited.clear();
						idTmp++;
						tmpSkips++;
					}
					// Seeds were linked, update the missing part of the long read
					if (linked != 0 && isLinkable) {
						if (correctedLR.empty()) {
							correctedLR = missingPart;
						} else {
							correctedLR = correctedLR + missingPart.substr(src.length());
						}
						src = tgt;
						posSrc = posTgt;
						if (seedsSkips > (seeds.size() - 1) - idSeed - 1) {
							seedsSkips = (seeds.size() - 1) - idSeed - 1;
						}
						nextReachable = tmpMissingPart.length() == 0 ? -1 : idTmp - 1;
						skippedSeeds = tmpMissingPart.length() == 0 ? 0 : tmpSkips - 1;
						firstSkippedSeed = skippedSeeds > 0 ? idSeed + 1 : -1;
					// Seeds couldn't be linked, skip a seed or split the corrected long read 
					} else {
						// Skip the target if the allowed number of skips isn't reached
						if (skippedSeeds < seedsSkips) {
							skippedSeeds++;
							if (firstSkippedSeed == -1) {
								firstSkippedSeed = idSeed;
							}
						// Split the corrected long read
						} else {
							if (!correctedLR.empty()) {
								if (posBeg > 0) {
									extendLeft(posBeg, correctedLR);
								}
								if ((int) (LRLen - posSrc - src.length()) > 0) {
									extendRight(LRLen - posSrc - src.length(), correctedLR);
								}
								fRes << ">" << LRId << "_" << fragments << endl << correctedLR << endl;
								fragments++;
							}
							
							// Update data
							if (firstSkippedSeed != -1) {
								idSeed = firstSkippedSeed;
							}
							curSeed = seeds[idSeed];
							src = curSeed.seq;
							posSrc = curSeed.pos;
							posBeg = posSrc;
							correctedLR = string();
							firstSkippedSeed = -1;
							skippedSeeds = 0;
							if (seedsSkips > (seeds.size() - 1) - idSeed - 1) {
								seedsSkips = (seeds.size() - 1) - idSeed - 1;
							}
						}
					}
					// Update the next seed to be processed
					idSeed = nextReachable != -1 ? nextReachable : idSeed + 1;
				}

				// Multiple seeds were mapped on the long read and were linked
				if (!correctedLR.empty()) {
					if (posBeg > 0) {
						extendLeft(posBeg, correctedLR);
					}
					if ((int) (LRLen - posSrc - src.length()) > 0) {
						extendRight(LRLen - posSrc - src.length(), correctedLR);
					}
					if (fragments == 1) {
						fRes << ">" << LRId << endl << correctedLR << endl;
					} else {
						fRes << ">" << LRId << "_" << fragments << endl << correctedLR << endl;
					}
				}

				// Only one seed was mapped on the long read
				if (seeds.size() < 2) {
					int srcLen = src.length();
					if (posBeg > 0) {
						extendLeft(posBeg, src);
					}
					if ((int) (LRLen - posSrc - srcLen) > 0) {
						extendRight(LRLen - posSrc - srcLen, src);
					}
					fRes << ">" << LRId << endl << src << endl;
				}
				
				outMtx.lock();
				cout << fRes.str();
				outMtx.unlock();
			}
			
			curLink++;
		}
	}
		
	void startCorrection(PgSAIndexStandard* index, unsigned maxorder, string tmpdir, unsigned seedsoverlap, unsigned minorder, unsigned maxbranches, unsigned maxseedsskips, unsigned nbThreads) {
		// Global variables
		pgsaIndex = index;
		minOrder = minorder;
		maxOrder = maxorder;
		maxBranches = maxbranches;
		seedsOverlap = seedsoverlap;
		maxSeedsSkips = maxseedsskips;
		tmpDir = tmpdir;
		
		// KMC database, for counting K-mers
		openDatabase(tmpDir + "/mers.db");
		
		// Prepare threads data
		vector<vector<string>> seeds;
		for (unsigned i = 0 ; i < nbThreads ; i++) {
			seeds.push_back(vector<string>());
		}
		ifstream f(tmpDir + "/seeds");
		string line;
		unsigned curT = 0;
		while(getline(f, line)) {
			seeds[curT % nbThreads].push_back(line);
			curT++;
		}
		
		// Launch threads
		vector<future<void>> threads;
		for (unsigned i = 0 ; i < nbThreads ; i++) {
			vector<string> tmpv = seeds[i];
			threads.push_back(async(launch::async, [tmpv]() mutable {
				generateCLRs(tmpv);
			}));
		}
		
		// Get threads results
		for (future<void> &t: threads) {
			t.get();
		}
		  
		closeDatabase();
	}
}
