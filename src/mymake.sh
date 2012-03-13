g++ -Wall -I/home/analysis/rsk/moore_anal/gsl/include/ -c multinomial.c
g++   -g -O3 -Wall -c model.c
g++  -I/home/analysis/rsk/moore_anal/gsl/include/ -I/home/linux/libs/sparsehash/1.11/include  -I/sw/compbio/samtools/build/samtools-0.1.16/  -g -O3 -Wall -c fragment_coverage.c
g++   -I/home/linux/libs/sparsehash/1.11/include  -I/sw/compbio/samtools/build/samtools-0.1.16/  -g -O3 -Wall -c bam_coverage.c

g++   -L/home/analysis/rsk/moore_anal/gsl/lib -lgsl -lgslcblas -lm -L/home/linux/libs/sparsehash/1.11/lib -I/home/linux/libs/sparsehash/1.11/include  -I/sw/compbio/samtools/build/samtools-0.1.16/ -g -O3 -Wall bam_coverage.o multinomial.o fragment_coverage.o model.o fgetopt.c qev.c -o qev -lz -L/sw/compbio/samtools/build/samtools-0.1.16/ -lbam


