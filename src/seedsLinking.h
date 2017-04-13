#ifndef SLRGEN_H
#define	SLRGEN_H

#include "test/testdata.h"
#include "index/PgSAIndexFactory.h"

using namespace PgSAReadsSet;
using namespace PgSAIndex;

namespace SLRgen {
    void generateSLR(PgSAIndexStandard* index, int KLEN, char* tolink, string tplName, int seedsoverlap, int minoverlap, int backtracks, int seedskips);
}

#endif	/* SLRGEN_H */

