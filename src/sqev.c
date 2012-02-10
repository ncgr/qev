#include <stdio.h>
#include <iostream>
#include <google/sparse_hash_map>
#include <string.h>
#define _cplusplus
#include "bam.h"
#include "sam.h"
#include "faidx.h"
#include <list>
#include <cstdio>
#include "fragment_coverage.h"
#include <fstream>

using namespace std;
using google::sparse_hash_map;

using __gnu_cxx::hash;

typedef struct {
        int beg, end;
        samfile_t *in;
} tmpstruct_t;

/*for the time being writing code in this file, ideally a scan stat library will be created*/








int main(int argc, char *argv[])
{
        tmpstruct_t tmp;
	bamFile fp;
	faidx_t *fai;
	if (argc != 3) {
		fprintf(stderr, "Usage: qev <in.bam> <in.fna>\n");
		return 1;
	}
        fp = bam_open(argv[1], "rb");
	fai = fai_load(argv[2]);	

	if (fp == 0) {
		fprintf(stderr, "Fail to open BAM file %s\n", argv[1]);
		return 1;
	}
        tmp.beg = 0; tmp.end = 0x7fffffff;
        tmp.in = samopen(argv[1], "rb", 0);
	ofstream dnull;
	dnull.open("/dev/null");

        if (tmp.in == 0) {
                fprintf(stderr, "Fail to open BAM file %s\n", argv[1]);
                return 1;
        } else {
                tmp.beg = 0; tmp.end = 0x7fffffff;

                bam_index_t *idx;
                idx = bam_index_load(argv[1]); // load BAM index
                if (idx == 0) {
                        fprintf(stderr, "BAM indexing file is not available.\n");
                        return 1;
                }
		bam_header_t *hin;
	        hin = bam_header_read(fp);
	
		for (int j = 0; j < hin->n_targets; ++j){
			hash_list_t *hash_list = new hash_list_t();
	                bam_fetch(tmp.in->x.bam, idx, j, tmp.beg, tmp.end, hash_list, fetch_func);
			remove_singlets(hash_list->plist);
			printf(">%s ", hin->target_name[j]);
			calculate_transcript(hash_list->plist,hash_list->bad_plist, tmp.in->header->target_len[j],&dnull);
	        	char *seq;
			int len;
        		seq = faidx_fetch_seq(fai, tmp.in->header->target_name[j], 0,  tmp.in->header->target_len[j] , &len);
			cout << seq << "\n";
			delete hash_list;
		}
	}
	return 0;
}

