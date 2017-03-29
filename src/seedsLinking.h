/*
 * File:   PgSATests.h
 * Author: Tomek
 *
 * Created on 6 luty 2014, 18:19
 */

#ifndef PGSATESTS_H
#define	PGSATESTS_H

#include "test/testdata.h"
#include "index/PgSAIndexFactory.h"

using namespace PgSAReadsSet;
using namespace PgSAIndex;

namespace PgSATest {

    void runTest(PgSAIndexStandard* index, int _repeat, int _testKmersCount, unsigned short k, bool pairScrambled, bool filterTTTs, bool byPosition, int KLEN, char* tolink, char* tplName, int seedsoverlap, int minoverlap, int backtracks, int seedskips);
}

#endif	/* PGSATESTS_H */

