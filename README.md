# HG-CoLoR
HG-CoLoR (Hybrid Graph for the error Correction of Long Reads) is a hybrid method for the
error correction of long reads that follows the main idea from NaS to produce corrected
long reads from assemblies of related accurate short reads.

HG-CoLoR however, instead of aligning all the short reads against each other, focuses on
a seed-and-extend approach based on a hybrid structure between a de Bruijn graph and an
overlap graph, built from the short reads. This hybrid graph allows to compute perfect
overlaps of variable length between the short reads' k-mers, and is used to extend and
link together the seeds, which are short reads that align correctly on the input long
reads, using them as anchor points. The corrected long reads are thus produced by directly
assembling the k-mers of the short reads during the graph traversal, without using any
other proper assembly tool.

Pre-requisites
--------------

  - A Linux based operating system.
  - Shell tool GNU Parallel available through your PATH environment variable (https://www.gnu.org/software/parallel/).
  - Emboss binaries accessible through your PATH environment variable (http://emboss.sourceforge.net/download/).
  - KMC3 binaries (kmc, kmc_tools and kmc_dump) accessible through your PATH environment variable (https://github.com/refresh-bio/KMC).
  - QuorUM binary accessible through your PATH environment variable (https://github.com/gmarcais/Quorum).
  - PgSA directory accessible somewhere on your computer (https://github.com/kowallus/PgSA).
  
Dependencies
--------------

The blasr binary comes from the blasr software. Copyright notice is given in the file
bin/blasr-license.
The PgSAgen.cpp file was copied, and the SLRgen.cpp, seedsLinking.cpp, and seedsLinking.h
files were adapted from the PgSA sources.
  
Installation
--------------

  1. Go to the PgSA directory, and compile the PgSA library:  
  `make build CONF=pgsalib`
  2. Go back to the HG-CoLoR directory and run:             
  `make PGSA_PATH=/absolute/path/to/your/PgSA/folder/`
  
Running HG-CoLoR
--------------

To run HG-CoLoR, run the following command:       
`./HG-CoLoR --longreads LR.fasta --shortreads SR.fastq --out result.fasta --tmpdir tmp_directory`

### Input

  - LR.fasta:       fasta file of long reads, one sequence per line.
  - SR.fastq:       fastq file of short reads.
    Warning: only one file must be provided.
    If using paired reads, please concatenate them into one single file.

  - result.fasta:   fasta file where to output the corrected long reads.
  - tmp_directory directory where to store the temporary files.

### Options

      --kmer:           k-mer size for the graph construction (default: 64).
      --solid:          Minimum number of occurrences to consider a k-mer as solid (default: 5).
                        This parameter should be raised accordingly to the short reads coverage and accuracy.
                        Its default value is adapted for a 50x coverage set of short reads with a 1% error rate.
      --seedsoverlap:   Minimum overlap length to allow the merging of two overlapping seeds (default: k-1).
      --minoverlap:     Minimum overlap length to allow the exploration of an edge of the graph (default: k-5).
      --backtracks:     Maximum number of backtracks (default: 1,000).
                        Raising this parameter will result in less fragmented corrected long reads.
                        However, it will also increase the runtime, and may create chimeric linkings between the seeds.
      --seedskips:      Maximum number of seed skips (default: 5).
      --bestn:          Top alignments to be reported by BLASR (default: 30).
                        This parameter should be raised accordingly to the short reads coverage.
                        Its default value is adapted for a 50x coverage of short reads.
      --kmcmem:         Maximum amount of RAM for KMC, in GB (default: 12)
      --nproc:          Number of processes to run in parallel (default: number of cores).
      --help:           Print a help message.

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

If you use HG-CoLoR in your research, please cite:

Pierre Morisse, Thierry Lecroq and Arnaud Lefebvre. HG-CoLoR: Hybrid-Graph
for the error Correction of Long Reads, Actes des Journées Ouvertes Biologie
Informatique et Mathématiques, 2017.

Contact
--------------

You can report problems and bugs to pierre[dot]morisse2[at]univ-rouen[dot]fr
