#include "seedsMerging.h"

int* computeBacktrackTable(std::string s) {
	int* T = (int*) malloc(s.length() * sizeof(int));
	int cnd = 0;
	T[0] = -1;
	T[1] = 0;
	int pos = 2;
	while (pos < s.length()) {
		if (s[pos - 1] == s[cnd]) {
			T[pos] = cnd + 1;
			pos += 1;
			cnd += 1;
		} else if (cnd > 0) {
			cnd = T[cnd];
		} else {
			T[pos] = 0;
			pos += 1;
		}
	}

	return T;
}

int overlapLength(std::string s1, std::string s2) {
	if (s1.length() > s2.length()) s1 = s1.substr(s1.length() - s2.length());

	int* T = computeBacktrackTable(s2);

	int m = 0;
	int i = 0;
	while (m + i < s1.length()) {
		if (s2[i] == s1[m + i]) {
			i += 1;
		} else {
			m += i - T[i];
			if (i > 0) i = T[i];
		}
	}
	
	free(T);
	return i;
	}

void mergeOverlappingPosSeeds(std::vector<seed_t> &seeds, unsigned minOverlap) {
  unsigned i = 0;
  unsigned j = 1;
  std::string s1;
  std::string s2;
  // Iterate through the seeds list
  while (i < seeds.size() - 1 && j < seeds.size()) {
	// Check if the two current seeds overlap
    if (seeds[i].pos != -1 && seeds[j].pos >= seeds[i].pos && seeds[j].pos <= seeds[i].pos + seeds[i].alen) {
      int b1 = seeds[j].pos - seeds[i].pos;
      
      s1 = seeds[i].seq.substr(b1);
      s2 = seeds[j].seq.substr(0, s1.length());
      // If the overlap is long enough, and if their overlapping sequences match, merge the two seeds
      if (seeds[j].pos + seeds[j].alen > seeds[i].pos + seeds[i].alen && s1.length() >= minOverlap && s1 == s2) {
        int cplen = seeds[j].pos + seeds[j].alen - seeds[i].pos - seeds[i].alen;
        seeds[i].seq = seeds[i].seq + seeds[j].seq.substr(seeds[j].alen - cplen);
        seeds[i].alen = seeds[i].alen + cplen;
        seeds[i].matches = seeds[i].matches + (seeds[j].matches / seeds[j].alen) * (cplen);
        seeds.erase(seeds.begin() + j);
      // Othewhise, keep the seed with the best alignment score
      } else if (seeds[i].matches < seeds[j].matches) {
        seeds.erase(seeds.begin() + j);
      } else {
        seeds.erase(seeds.begin() + i);
      }
    // The two seeds don't overlap, move to the next ones
    } else {
		i = j;
		j++;
	}
  }
}

void mergeOverlappingSeqSeeds(std::vector<seed_t> &seeds, int minOverlap) {
  unsigned i = 0;
  unsigned j = 1;
  // Iterate through the seeds list
  while (i < seeds.size() - 1 && j < seeds.size()) {
	int overlap = overlapLength(seeds[i].seq, seeds[j].seq);
	// If the two current seeds overlap on a sufficient length, merge them
	if (overlap >= minOverlap) {
		int cplen = seeds[j].alen - overlap;
        seeds[i].seq = seeds[i].seq + seeds[j].seq.substr(overlap);
		seeds[i].alen = seeds[i].alen + cplen;
		seeds[i].matches = seeds[i].matches + seeds[j].matches / 2;
		seeds.erase(seeds.begin() + j);
	} else {
		i = j;
		j++;
	}
  }
}

int getTemplateLength(char* tpl, int len) {
	int i = len - 1;
	while (i >= 0 && tpl[i] != '_') {
		i--;
	}
	i++;
	return atoi(tpl + i);
}

std::vector<seed_t> readAlignmentFile(char* alFile) {
	FILE* f = fopen(alFile, "r");
	if (f == NULL) {
		perror("");
		exit(1);
	}
	
	int posT;
	int rlen;
	int tlen = -1;
	int matches;
	std::string seq;
	char line[4096];
	char* t;
	std::vector<seed_t> seeds;
  
	// Iterate through the lines of the file, and get the interesting fields
	while (fgets(line, sizeof(line), f) != NULL) {
		t = strtok(line, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		if (tlen == -1) {
			tlen = getTemplateLength(t, strlen(t));
		}
		t = strtok(NULL, "\t");
		posT = atoi(t);
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		seq = std::string(t);
		rlen = seq.length();
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		matches = atoi(t + 5);
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		
		seeds.push_back({posT, rlen, tlen, matches, seq});
	}
	
	fclose(f);
	
	return seeds;
}

std::vector<seed_t> processSeeds(char* alFile, unsigned minOverlap) {
	std::vector<seed_t> seeds = readAlignmentFile(alFile);
  
	// Sort the seeds according to their mapping positions
	std::sort(seeds.begin(), seeds.end());
	// Merge the seeds
	mergeOverlappingPosSeeds(seeds, minOverlap);
	mergeOverlappingSeqSeeds(seeds, minOverlap);

	return seeds;
} 
