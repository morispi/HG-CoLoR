# HG-CoLoR
HG-CoLoR (Hybrid method based on a variable-order de bruijn Graph for the error Correction of Long Reads)
is a hybrid method for the error correction of long reads that both aligns the short reads to the long reads,
and uses a variable-order de Bruijn graph, in a seed-and-extend approach. The seeds, found by aligning
the short reads to the long reads, are used as anchor points on the variable-order de Bruijn graph, built
from the short reads, which is traversed in order to find paths allowing to link seeds together. Such paths
between seeds dictate corrections for the missing part of the long reads, that are not covered by seeds.

Requirments
--------------

  - A Linux based operating system.
  - Python3.
  - Emboss binaries accessible through your PATH environment variable (http://emboss.sourceforge.net/download/).
  - KMC3 binaries accessible through your PATH environment variable (https://github.com/refresh-bio/KMC).
  - QuorUM binary accessible through your PATH environment variable (https://github.com/gmarcais/Quorum).
  
Dependencies
--------------

The blasr binary comes from the blasr software. Copyright notice is given in the file
bin/blasr-license.
  
Installation
--------------
  ```bash
  git clone https://github.com/pierre-morisse/HG-CoLoR
  git submodule init
  git submodule update
  cd KMC/ && make -j
  cd ../PgSA/ && make build CONF=pgsalib
  cd .. && make
  ```
  
Running HG-CoLoR
--------------

To run HG-CoLoR, run the following command:

`./HG-CoLoR --longreads LR.fasta --shortreads SR.fastq --out resultPrefix -K maxK`

### Input

  - LR.fasta:       fasta file of long reads, one sequence per line.
  - SR.fastq:       fastq file of short reads.
    Warning: only one file must be provided.
    If using paired reads, please concatenate them into one single file.

  - resultPrefix:   Prefix of the fasta files where to output the corrected, trim and split long reads.
  - maxK:           Maximum K-mer size of the variable-order de Bruijn graph.

### Output format

The corrected reads are output in fasta format, with one sequence per line. The header of each corrected read
consists of 5 components, as follows:

`>id_len_seedsBases_graphBases_rawBases`

where

  - `id` is the original read header
  - `len` is the original read length
  - `seedsBases` is the number of bases of the corrected long read coming from seeds
  - `graphBases` is the number of bases of the corrected long read coming from the traversals of the variable-order de Bruijn graph
  - `rawBases` is the number of (uncorrected) bases of the corrected long read, coming from the original, raw long read

### Options

      --minorder INT, -k INT:       Minimum order of the variable-order de Bruijn graph (default: K/2).
      --solid INT, -S INT:          Minimum number of occurrences to consider a k-mer as solid (default: 1).
                                    This parameter should be set accordingly to the short reads coverage and accuracy,
                                    and to the chosen maximum order of the graph.
                                    It should only be increased when using high coverage of short reads, or a small maximum order.
      --seedsoverlap INT, -o INT:   Minimum overlap length to allow the merging of two overlapping seeds (default: maxorder - 1).
      --seedsdistance INT, -d INT:   Maximum distance to consider two consecutive seeds for merging (default: 10).
      --branches INT, -b INT:       Maximum number of branches exploration (default: 1,500).
                                    Raising this parameter will result in less split corrected long reads.
                                    However, it will also increase the runtime, and may create chimeric links between the seeds.
      --seedskips INT, -s INT:      Maximum number of seed skips (default: 5).
      --mismatches INT, -m INT:     Allowed mismatches when attempting to link two seeds together (default: 3).
      --bestn INT, -n INT:          Top alignments to be reported by BLASR (default: 50).
                                    This parameter should be set accordingly to the short reads coverage.
                                    Its default value is adapted for a 50x coverage of short reads.
                                    It should be decreased with higher coverage, and increased with lower coverage.
      --nproc INT, -j INT:          Number of processes to run in parallel (default: number of cores).
      --tmpdir STRING, -t STRING:   Path where to store the directory containing temporary files (default: working directory)
      --kmcmem INT, -r INT:         Maximum amount of RAM for KMC, in GB (default: 12)
      --help, -h:                   Print this help message.

### Short reads coverage and accuracy

HG-CoLoR default parameters are adapted for a 50x coverage set of short reads with a 1% error rate. Please modify the parameters, in particular the --solid and --bestn ones,
as indicated above if using a set of short reads with a much higher coverage and/or a highly different error rate.
      
Notes
--------------

HG-CoLoR has been developed and tested on x86-64 GNU/Linux.          
Support for any other platform has not been tested.

Authors
--------------

Pierre Morisse, Thierry Lecroq and Arnaud Lefebvre.

Reference
--------------

The HG-CoLoR paper is accepted for publication in Bioinformatics under the title "Hybrid correction of highly noisy long reads using a variable-order de Bruijn graph"!

Contact
--------------

You can report problems and bugs to pierre[dot]morisse2[at]univ-rouen[dot]fr
