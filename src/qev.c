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
#include "fgetopt.h"
#include "bam_coverage.h"

using namespace std;
using google::sparse_hash_map;

using __gnu_cxx::hash;

typedef struct {
        int beg, end;
        samfile_t *in;
} tmpstruct_t;

int main(int argc, char *argv[])
{
        tmpstruct_t tmp;
	tmp.in=0;
	bamFile fp = 0;
	faidx_t *fai = 0;
	 bam_index_t *idx=0;

   static struct option long_options[] =
        {
                        {"help",0,0,'h'},  //help 
                        {"bam", 1, 0, 'b'}, // the bam file for input required
                        {"fasta", 1, 0, 'f'},// the fasta file for input required
                        {"transcript", 1, 0, 't'}, // optional only get one transcript
                        {"log", 1, 0, 'l'}, // optional output log of the coverage
			{0,0,0,'h'}       //print help if no commands ****this element of the array must be at the end.
        };
    char usageMessage[] =
    "\n\nOptions for this program include\n\
    -h, --help,           Print this message.\n\
    -b, --bam,           required bam indexed bam file for input.\n\
\n\
\n\
\n\
    Optional arguments \n\
    -f, --fasta    generates output as a fasta, with a fasta file for sequence input \n\
    -t, --transcript evaluates the quality for a single transcript \n\
    -l, --log outputs additional log information \n\n";

    int option_index=0;
    int optflag=0;
    char *transcript=0;
    char *logfilename=0;


   

    while(  (optflag = fgetopt_long ( argc, argv , "b:f:t:l:h" , long_options , &option_index )) !=EOF ){
       printf("Option %c Argument %s\n",optflag,optarg);
	printf ("option %s\n", long_options[option_index].name);
        switch(optflag){
        case 'b':
                cerr << "Bam file -b " << optarg << "\n";
		fp = bam_open(optarg, "rb");
		tmp.in = samopen(optarg, "rb", 0);
                idx = bam_index_load(optarg); // load BAM index

		if (fp == 0 || tmp.in ==0) {
        	        fprintf(stderr, "Fail to open BAM file or index of BAM file%s\n", optarg);
	                return 1;
		}
        	
		break;
        case 'f':
                cerr << "fasta file -f " << optarg << "\n";
		fai = fai_load(optarg);
		if (fai == 0) {
                        fprintf(stderr, "Fail to open fasta file %s\n", optarg);
                        return 1;
                }
		break;
        case 't':
                cerr << "transcript name -t " << optarg << "\n";
		transcript = new char[strlen(optarg)+1];
		strcpy(transcript,optarg);
		break;
        case 'l':
                cerr << "log file -l " << optarg << "\n";
		logfilename = new char[strlen(optarg)+1];
                strcpy(logfilename,optarg);

		break;
	default:
	 case ':':
        case '?':
        case 'h':
                printf(usageMessage);
                return 1;
	}
	}

	if(!fp){
		printf(usageMessage);
                return 1;
	}
	ofstream logfile;
        if(!logfilename){
                logfile.open("/dev/null");
        }else{
                logfile.open(logfilename);
        }



        tmp.beg = 0; tmp.end = 0x7fffffff;

	bam_header_t *hin;
	hin = bam_header_read(fp);
	
	for (int j = 0; j < hin->n_targets; ++j){
		if( !transcript || strcmp(hin->target_name[j], transcript)==0){
			if(fai){
				printf(">");
			}
			printf("%s ", hin->target_name[j]);
			hash_list_t *hash_list = new hash_list_t();
		        bam_fetch(tmp.in->x.bam, idx, j, tmp.beg, tmp.end, hash_list, fetch_func);
			remove_singlets(hash_list->plist);
			logfile << ">" << hin->target_name[j] << "\n";
			calculate_transcript_shape_coverage(hash_list->plist,hash_list->bad_plist, tmp.in->header->target_len[j],&logfile);
			
	        	char *seq;
			int len;
			if(fai){
	      			seq = faidx_fetch_seq(fai, tmp.in->header->target_name[j], 0,  tmp.in->header->target_len[j] , &len);
				cout << seq << "\n";
			}
			delete hash_list;
		}
	}
	return 0;
}

