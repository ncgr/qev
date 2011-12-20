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
#include <gsl/gsl_cdf.h>
#include <math.h>
#include "multinomial.h"
#include "fragment_coverage.h"

using namespace std;
using google::sparse_hash_map;

using __gnu_cxx::hash;
class model_t{
	public:
	double *p;
        double *E;
        double *pI;
        unsigned int *EI;
	model_t(unsigned long size, unsigned long *insert_counts, unsigned long counts);
	~model_t();	
};


/*
creates a model probability model for testing
*/
model_t::model_t(unsigned long size, unsigned long *insert_counts, unsigned long counts){
        p = new double[size];
        E = new double[size];
        pI = new double[size];
        EI = new unsigned int[size];
        memset(p,0,sizeof(double)*size);
		
	for (unsigned int i = 0;i<= size; i++){
		if(insert_counts[i]>0){
			 unsigned int FL = i;//fragment length
			double temp_p=0.0;	
			if(size < 2*FL){
			//model I: inserts will always overlap
		                for (unsigned int inc = 1; inc <= size; inc ++){
		                        if(inc < size - FL){ //if size==4, FL==3
                		                temp_p=((double)inc)/((double)(size-FL+1));
		                        }else if(inc <= FL){
                		                temp_p = 1.0 ;
		                        }else{  //
						temp_p=((double)(size-inc+1))/((double)(size-FL+1));
					}
                		        E[inc-1]+=((double)insert_counts[i])*temp_p;
		                }
			}else {
                                for (unsigned int inc = 1; inc <= size; inc ++){
                                        if(inc < FL){     
                                                 temp_p=((double)inc)/(double)(size-FL+1);
                                        }else if(inc <= size - FL ){
                                                temp_p = ((double)FL)/(double)(size-FL+1);
                                        }else{ 
                                                temp_p=((double)(size-inc+1))/(double)(size-FL+1);
                                        }
                                        E[inc-1]+=((double)insert_counts[i])*temp_p;
                                }

			}
		}
	}

	if(counts>0){
		for (unsigned int i = 0;i< size; i++){
			p[i]=E[i]/(double)counts;
			cout << "model_E\t" << i<< "\t" << E[i] << "\n";
		}
	}else{
		for (unsigned int i = 0;i< size; i++){
                        p[i]=0.0;
                        E[i]=0.0;
                }
	}
}
model_t::~model_t(){
                delete p;
                delete E;
                delete pI;
                delete EI;
}


void calculate_transcript(list<pair_t *> *plist, list<pair_t *> *bad_plist, unsigned int target_len){
        //process contig calculating paired insert coverage.
        unsigned long *pcov=new unsigned long[target_len];
        unsigned long *fragment_hist=new unsigned long[target_len+1];
        memset(pcov,0,sizeof(unsigned long)*target_len);
        memset(fragment_hist,0,sizeof(unsigned long)*(target_len+1));
        //Clean up hash and list for next set of pairs.
        double mean_insert=0.0;
        unsigned long counts=0;
        //calculate coverage, and mean insert length.
        calculate_fragment_coverage(plist,pcov,&mean_insert,&counts,fragment_hist,target_len);
        model_t *model = new model_t(target_len, fragment_hist, counts);
        //calculate the minimum probability of any position in the transcript.
        double min_binom_p=1.0;
        unsigned int min_ink=0;
        double p_total=0.0;
        for(unsigned int ink=0;ink<target_len;++ink){
        	double P=0.0;
                if(model->p[ink] < 0.0 || model->p[ink] > 1.0){
                	cout << "p out-of-bounds " << model->p[ink] << " at " << ink << "\n";
                        P=-1;
                }else{
	                P = gsl_cdf_binomial_P(pcov[ink],model->p[ink],counts);
                }
        	if(min_binom_p>P){
                	min_binom_p = P;
	                min_ink =  ink;
                }
                p_total+=model->p[ink];
        }

        //set the relative lower bound for subsequently calculating the confidence in the transcript.
        unsigned int N=0;
        int bad_count=0;
        for(unsigned int ink=0;ink<target_len;++ink){
        	N+=pcov[ink];
                model->pI[ink]=model->p[ink]/p_total;
                //set relative lower bound
                model->EI[ink]=(unsigned int) (model->E[ink]*(pcov[min_ink]/model->E[min_ink])+.5);
	        if( model->p[ink]<0.0 || model->p[ink]>1.0){
        	        bad_count++;
                }
		cout << "pcov_P\t" << ink << "\t" << pcov[ink] <<"\n";
	}
        double res = 0;
        if(N>0 || bad_count==0){
        	res = cdf_multinomial_P(target_len,N,model->pI,model->EI);
        }else{
                res = -1;
        }
        cout << "len " << target_len <<" frag_c "<< counts << " tot_Nucs " << N << " ave_frag_len " << mean_insert << " pos " << min_ink <<" minPcdf " << min_binom_p << " model_p " << model->p[min_ink] << " model_expec " << model->E[min_ink] << " m_cdf " << res << " pcov " << pcov[min_ink] << " prop " << plist->size() << " imp " << bad_plist->size() << " p_tot " << p_total <<"\n";
        //cout << "\n";
        delete plist;
	delete bad_plist;
        delete pcov;
        delete fragment_hist;
        delete model;
}
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

