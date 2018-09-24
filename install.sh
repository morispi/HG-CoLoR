#!/bin/bash

git submodule init
git submodule update
cd KMC/
make -j
cd ../PgSA/
make build CONF=pgsalib
make build CONF=pgsagen
cd ..
make
