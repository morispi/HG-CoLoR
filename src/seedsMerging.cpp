#include "seedsMerging.h"

char* strToUpper(char* str, int len) {
	for (int i = 0; i < len; i++) {
		str[i] = toupper(str[i]);
	}
	return str;
}

int* computeBacktrackTable(char* s) {
	int* T = (int*) malloc(strlen(s) * sizeof(int));
	int cnd = 0;
	T[0] = -1;
	T[1] = 0;
	unsigned pos = 2;
	while (pos < strlen(s)) {
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

int overlapLength(char* s1, char* s2) {
	char tmp[strlen(s2) + 1];
	if (strlen(s1) > strlen(s2)) {
		 strncpy(tmp, s1 + strlen(s1) - strlen(s2), strlen(s2) + 1);
	} else {
		strncpy(tmp, s1, strlen(s1) + 1);
	}

	int* T = computeBacktrackTable(s2);

	unsigned m = 0;
	unsigned i = 0;
	while (m + i < strlen(tmp)) {
		if (s2[i] == tmp[m + i]) {
			i += 1;
		} else {
			m += i - T[i];
			if (i > 0) i = T[i];
		}
	}

	return i;
}

int compareAlignmentsPos(const void* a1, const void* a2) {
  seed_t* t1 = (seed_t*) a1;
  seed_t* t2 = (seed_t*) a2;
  return t1->pos - t2->pos;
}

seed_t* mergeOverlappingPosSeeds(int* nb, seed_t* seeds, unsigned minOverlap) {
  int n = 0;
  int i = 0;
  int j = 1;
  char* s1 = NULL;
  char* s2 = NULL;
  // Parcours de la liste, calcul des chevauchement, et élimination des reads avec un plus faible nombre de matches
  while (i < *nb - 1 && j < *nb) {
    if (seeds[i].pos != -1 && seeds[j].pos >= seeds[i].pos && seeds[j].pos <= seeds[i].pos + seeds[i].alen) {
      int b1 = seeds[j].pos - seeds[i].pos;
      s1 = (char*) malloc(seeds[i].alen - b1 + 1);
      s2 = (char*) malloc(seeds[i].alen - b1 + 1);
      strncpy(s1, seeds[i].seq + b1, seeds[i].alen - b1 + 1);
      strncpy(s2, seeds[j].seq, seeds[i].alen - b1);
      s2[seeds[i].alen - b1] = '\0';
      if (seeds[j].pos + seeds[j].alen > seeds[i].pos + seeds[i].alen && strlen(s1) >= minOverlap && strcmp(s1, s2) == 0) {
        // fusionner
        int cplen = seeds[j].pos + seeds[j].alen - seeds[i].pos - seeds[i].alen;
        char* tmp = (char*) malloc(seeds[i].alen + cplen + 1);
        strncpy(tmp, seeds[i].seq, seeds[i].alen + 1);
        strncat(tmp, seeds[j].seq + seeds[j].alen - cplen, cplen);
        seeds[i].alen = seeds[i].alen + cplen;
        seeds[i].matches = seeds[i].matches + (seeds[j].matches / seeds[j].alen) * (cplen);
        free(seeds[i].seq);
        seeds[i].seq = tmp;
        free(seeds[j].seq);
        seeds[j].seq = NULL;
        seeds[j].pos = -1;
        j++;
      } else if (seeds[i].matches < seeds[j].matches) {
        seeds[j].pos = -1;
        free(seeds[j].seq);
        seeds[j].seq = NULL;
        j++;
      } else {
        seeds[i].pos = -1;
        free(seeds[i].seq);
        seeds[i].seq = NULL;
        i = j;
        j++;
      }
      n++;
      free(s1);
      s1 = NULL;
      free(s2);
      s2 = NULL;
    } else {
		i = j;
		j++;
	}
  }
  
  // Pas de chevauchement, retour de la liste initiale.
  if (n == 0) {
    return seeds;
  }
  
  // Création et initialisation de la nouvelle liste, mise à jour de *nb,
  // et retour.
  seed_t* res = (seed_t*) malloc((*nb - n) * (sizeof(*res)));
  i = 0;
  j = 0;
  while (i < *nb) {
    if (seeds[i].pos != -1) {
      res[j] = seeds[i];
      j++;
    }
    i++;
  }
  *nb = *nb - n;
  return res;
}

seed_t* mergeOverlappingSeqSeeds(int* nb, seed_t* seeds, int minOverlap) {
  int n = 0;
  int i = 0;
  int j = 1;
  // Parcours de la liste, calcul des chevauchement, et élimination des reads avec un plus faible nombre de matches
  while (i < *nb - 1 && j < *nb) {
	int overlap = overlapLength(seeds[i].seq, seeds[j].seq);
	if (overlap >= minOverlap) {
		int cplen = seeds[j].alen - overlap;
		char* tmp = (char*) malloc(seeds[i].alen + cplen + 1);
		strncpy(tmp, seeds[i].seq, seeds[i].alen + 1);
		strncat(tmp, seeds[j].seq + overlap, cplen);
		seeds[i].alen = seeds[i].alen + cplen;
		seeds[i].matches = seeds[i].matches + seeds[j].matches / 2;
		free(seeds[i].seq);
		seeds[i].seq = tmp;
		free(seeds[j].seq);
		seeds[j].seq = NULL;
		seeds[j].pos = -1;
		n++;
		j++;
	} else {
		i = j;
		j++;
	}
  }
  
  // Pas de chevauchement, retour de la liste initiale.
  if (n == 0) {
    return seeds;
  }
  
  // Création et initialisation de la nouvelle liste, mise à jour de *nb,
  // et retour.
  seed_t* res = (seed_t*) malloc((*nb - n) * (sizeof(*res)));
  i = 0;
  j = 0;
  while (i < *nb) {
    if (seeds[i].pos != -1) {
      res[j] = seeds[i];
      j++;
    }
    i++;
  }
  *nb = *nb - n;
  return res;
}

int countAlignments(char* file) {
	FILE*f = fopen(file, "r");
	
	if (f == NULL) {
		perror("");
		exit(1);
	}
	
	char line[4096];
	int matches;
	
	while (fgets(line, sizeof(line), f) != NULL) {
		matches++;
	}
	
	fclose(f);
	return matches;
}

int getTemplateLength(char* tpl, int len) {
	int i = len - 1;
	while (i >= 0 && tpl[i] != '_') {
		i--;
	}
	i++;
	return atoi(tpl + i);
}

void readAlignmentFile(char* file, int n_seeds, seed_t** seeds) {
	FILE* f = fopen(file, "r");
	if (f == NULL) {
		perror("");
		exit(1);
	}
  
	char line[4096];
	char* t;
	seed_t* t_seeds = (seed_t*) malloc(n_seeds * sizeof(seed_t));
  
	// Indices pour le parcours des listes
	int i_seeds = 0;
  
	// Type, position et longueur de l'alignement, et longueur du read
	int posT;
	int rlen;
	int tlen = -1;
	int matches;
	char* seq;
  
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
		seq = (char*) malloc(strlen(t) + 1);
		strncpy(seq, t, strlen(t));
		seq[strlen(t)] = '\0';
		rlen = strlen(seq);
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		matches = atoi(t + 5);
		strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		t = strtok(NULL, "\t");
		
		t_seeds[i_seeds].pos = posT;
		t_seeds[i_seeds].alen = rlen;
		t_seeds[i_seeds].tlen = tlen;
		t_seeds[i_seeds].seq = seq;
		t_seeds[i_seeds].matches = matches;
		i_seeds++;
	}
	
	*seeds = t_seeds;
}

void freeAlignments(int size, seed_t* seeds) {
	for (int i = 0; i < size; i++) {
		free(seeds[i].seq);
	}
	free(seeds);
}

std::pair<int, seed_t*> processSeeds(char* tpl, unsigned minOverlap) {
	int n_seeds = countAlignments(tpl);
	seed_t* seeds;
  
	// Récupération des informations d'alignement du fichier
	readAlignmentFile(tpl, n_seeds, &seeds);
  
	// Tri des listes
	qsort(seeds, n_seeds, sizeof(seed_t), compareAlignmentsPos);
  
	// Fusion de seeds
	seeds = mergeOverlappingPosSeeds(&n_seeds, seeds, minOverlap);
	seeds = mergeOverlappingSeqSeeds(&n_seeds, seeds, minOverlap);

	return std::make_pair(n_seeds, seeds);
} 
