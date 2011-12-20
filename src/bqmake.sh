g++  -I/home/linux/libs/sparsehash/1.11/include  -I/sw/compbio/samtools/build/samtools-0.1.16/ -g -O2 -Wall -c fragment_coverage.c

g++  -lm -L/home/linux/libs/sparsehash/1.11/lib -I/home/linux/libs/sparsehash/1.11/include  -I/sw/compbio/samtools/build/samtools-0.1.16/ -g -O2 -Wall fragment_coverage.o bqev.c -o bqev -lz -L/sw/compbio/samtools/build/samtools-0.1.16/ -lbam


