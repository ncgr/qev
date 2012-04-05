#include <stdio.h>
#include <iostream>
#include <google/sparse_hash_map>
#include <string.h>
#define _cplusplus
#include "bam.h"
#include <list>
#include <cstdio>
#include <math.h>
#include "fragment_coverage.h"
#include <fstream>

using namespace std;
using google::sparse_hash_map;

using __gnu_cxx::hash;


#define HIST_BUFFER 1000000

void calculate_fragment_coverage(list<pair_t *> *plist, unsigned long *pcov, double *mean_insert, unsigned long * counts, unsigned int length){

        pair_t *pair=NULL;

        //process contig calculating paired insert coverage.
        memset(pcov,0,sizeof(long)*length);
        *mean_insert=0.0;
        *counts=0;
        //calculate coverage, and mean insert length.
        for(  list<pair_t *>::iterator i=plist->begin(); i != plist->end(); ++i){
                pair = *i;
                if(pair->first && pair->second){
                        unsigned int start, stop;
                        calc_span(&start,&stop,pair);
                        //alignments can go past the end of the read.
                        start = (start <= 0)? 0 : start;
                        stop = (stop >= length)? length-1 : stop;
                        for(unsigned int count=start; count<=stop;count++){

                                pcov[count]++;
                        }
                        //cout << stop - start + 1 << "\t" << start << "\t" << stop << "\n";
                        (*mean_insert)+=(double)(stop - start + 1);
                        (*counts)++;
                }
        }
        if((*counts)>0.0){
                (*mean_insert)/=(double)(*counts);
        }

}



int length_hist(char *fname, unsigned long *hist, double **mean_mean_coverage, unsigned long **fragment_length_hists, double *mean_pair_length, list<long *> *properVsplit_list){
	phash_t *pairs;
	phash_t *bad_pairs;
        list<pair_t *> *plist;
	list<pair_t *> *bad_plist;
        hash_list_t *hash_list;
	 hash_list = new hash_list_t();
        bamFile fp;
        fp = bam_open(fname, "rb");

        if (fp == 0) {
                fprintf(stderr, "Fail to open BAM file %s\n", fname);
                return 1;
        }
        //read target names.
        bam_header_t *hin;
        hin = bam_header_read(fp);

        bam_index_t *idx;
        idx = bam_index_load(fname); // load BAM index
        if (idx == 0) {
                fprintf(stderr, "BAM indexing file is not available.\n");
                return 1;
        }
	
        for (int j = 0; j < hin->n_targets; ++j){
		//These only need to be allocated once.
		if( mean_mean_coverage[hin->target_len[j]]==0){
			mean_mean_coverage[hin->target_len[j]] = new double[hin->target_len[j]];			
                	memset(mean_mean_coverage[hin->target_len[j]],0,sizeof(double)*hin->target_len[j]);
			mean_pair_length[hin->target_len[j]]=0.0;
			fragment_length_hists[hin->target_len[j]] = new unsigned long[hin->target_len[j]+1];
			memset(fragment_length_hists[hin->target_len[j]],0,sizeof(double)*hin->target_len[j]+1);

		}
		unsigned long *pcov = new unsigned long[hin->target_len[j]];
		plist = new list<pair_t *>();
		bad_plist = new list<pair_t *>();
		
                pairs = new phash_t();
		bad_pairs = new phash_t();
                hash_list->plist=plist;
                hash_list->hash=pairs;
		hash_list->bad_plist=bad_plist;
                hash_list->bad_hash=bad_pairs;
		

                hash_list->plist->clear();
                hash_list->hash->clear();
                hash_list->bad_plist->clear();
                hash_list->bad_hash->clear();



                bam_fetch(fp, idx, j, 0, 0x7fffffff, hash_list, fetch_func);
                remove_singlets(plist);
		remove_singlets(bad_plist);
		long *propVsplit=new long[2];
		
		propVsplit[0]=plist->size();
		propVsplit[1]=bad_plist->size();
		(*properVsplit_list).push_back(propVsplit);
		double mean_insert;
		unsigned long counts=0;
		calculate_fragment_coverage(plist,pcov,&mean_insert,&counts,hin->target_len[j]);
//		calculate_paired_unpaired_broken_pair();		

		calculate_fragment_hist(plist,fragment_length_hists[hin->target_len[j]], hin->target_len[j]);

//		cout << "counts "<< counts <<"\n";
		if( counts>0){
			mean_pair_length[hin->target_len[j]]+=mean_insert;

			for(unsigned long si=0;si<hin->target_len[j];si++){
				if(pcov[si]>0){
					mean_mean_coverage[hin->target_len[j]][si]+=((double)pcov[si])/(double)counts;
				}
			}
			 hist[hin->target_len[j]]++;
		}
		delete pcov;
					
	}
	return 0;
}

int main(int argc, char *argv[])
{
        if (argc < 3 ) {
                fprintf(stderr, "Usage: bqev <out> <*.bam>\n");
                return 1;
        }
	unsigned long *hist;
	double **mean_mean_coverage = new double* [HIST_BUFFER];
	memset(mean_mean_coverage,0,sizeof(double *)*HIST_BUFFER);
	unsigned long **fragment_hists = new unsigned long* [HIST_BUFFER];
	memset(fragment_hists,0,sizeof(unsigned long *)*HIST_BUFFER);
	double *mean_pair_length = new double [HIST_BUFFER];



        hist=(unsigned long *)malloc(sizeof(unsigned long)*HIST_BUFFER);
        memset(hist,0,sizeof(unsigned long)*HIST_BUFFER);
	list<long *> *properVsplit_list = new list<long *>;
		
        for(int i=2;i<argc;i++){
                cerr << argv[i] << "\n";
                length_hist(argv[i],hist,mean_mean_coverage, fragment_hists, mean_pair_length,properVsplit_list);
        }
	ofstream pVsf;
        pVsf.open(argv[1]);
	for(  list<long *>::iterator i=properVsplit_list->begin(); i != properVsplit_list->end(); ++i){
		long *t = *i;
		pVsf << t[0] << "\t" << t[1] <<"\n";
	}
	pVsf.close();
	for(int i=0;i<HIST_BUFFER;i++){
		if(hist[i]>0){
			double mpl=mean_pair_length[i]/(double)hist[i];
			cout << i << "\tcount\t" << hist[i] << "\tmpl\t" << mpl <<"\t" << mean_pair_length[i]<< "\n";
			if(hist[i]>0){
				for(int j=0;j<i;j++){
					double mean_mean;
					mean_mean =  mean_mean_coverage[i][j]/(double)hist[i];
					cout <<"mean_pcov\t" << j << " " << mean_mean << "\t";
					unsigned long fragment_counts=fragment_hists[i][j];
					cout  << "frag\t" << i  << "\t" << j  << "\t" << fragment_counts << "\n";
				}
			}
		}
	}
        return 0;
}

