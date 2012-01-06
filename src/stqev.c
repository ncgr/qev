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

using namespace std;
typedef struct {
        int beg, end;
        samfile_t *in;
} tmpstruct_t;

int main(int argc, char *argv[])
{
        tmpstruct_t tmp;
	bamFile fp;
	faidx_t *fai;
	if (argc != 4) {
		fprintf(stderr, "Usage: stqev <in.bam> <in.fna> <transcript_name>\n");
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

        if (tmp.in == 0) {
                fprintf(stderr, "Fail to open BAM file %s\n", argv[1]);
                return 1;
        } else {
                tmp.beg = 0; tmp.end = 0x7fffffff;

                int ref;
                bam_index_t *idx;
                idx = bam_index_load(argv[1]); // load BAM index
                if (idx == 0) {
                        fprintf(stderr, "BAM indexing file is not available.\n");
                        return 1;
                }

                bam_parse_region(tmp.in->header, argv[3], &ref,
                                 &tmp.beg, &tmp.end); // parse the region
                if (ref < 0) {
                        fprintf(stderr, "Invalid region %s\n", argv[3]);
                        return 1;
                }
		 hash_list_t *hash_list = new hash_list_t();
                bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, hash_list, fetch_func);
		remove_singlets(hash_list->plist);
		calculate_transcript(hash_list->plist,hash_list->bad_plist, tmp.in->header->target_len[ref]);
        	char *seq;
		int len;
        	seq = faidx_fetch_seq(fai, tmp.in->header->target_name[ref], 0,  tmp.in->header->target_len[ref] , &len);
		cout << seq << "\n";
	}
	return 0;
}

