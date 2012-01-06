#include <stdio.h>
#include <iostream>
#include <google/sparse_hash_map>
#include <string.h>
#define _cplusplus
#include "bam.h"
#include <list>
#include <cstdio>
#include "fragment_coverage.h"
#include <gsl/gsl_cdf.h>
#include <math.h>
#include "multinomial.h"
#include <ostream>

hash_list_t::hash_list_t(){
        plist = new list<pair_t *>();
        bad_plist = new list<pair_t *>();
        hash = new phash_t();
        bad_hash = new phash_t();
}

hash_list_t::~hash_list_t(){
        pair_t *pair;
        for(  phash_t::iterator i=hash->begin(); i != hash->end(); ++i){
                pair = i->second;
                bam_destroy1(pair->first);
                bam_destroy1(pair->second);
                free(pair);
        }
       for(  phash_t::iterator i=bad_hash->begin(); i != bad_hash->end(); ++i){
                pair = i->second;
                bam_destroy1(pair->first);
                bam_destroy1(pair->second);
                free(pair);
        }
	delete plist;
	delete bad_plist;
	delete hash;
	delete bad_hash;
}

void calc_span(unsigned int *start, unsigned int *stop, pair_t *pair){
	unsigned int x[4];
	x[0] = pair->first->core.pos;
	x[1] = pair->second->core.pos;
	x[2] =bam_calend(&pair->first->core,  bam1_cigar(pair->first)); 
	x[3] = bam_calend(&pair->second->core,  bam1_cigar(pair->second));
	
	*start=x[0];
	*stop=x[1];
	for(int i=0; i<4;i++){
		if(*start>x[i]){
			*start=x[i];
		}
		if(*stop<x[i]){
			*stop=x[i];
		}
	}
}

int fetch_func( const bam1_t *b, void *data)
{
	uint32_t end;
	uint32_t end2;
	hash_list_t *hash_list = (hash_list_t *)data;
	phash_t *buf = hash_list->hash;
        list<pair_t *> *plist = hash_list->plist;
        phash_t *bad_buf = hash_list->bad_hash;
        list<pair_t *> *bad_plist = hash_list->bad_plist;


	pair_t *pair=NULL;
	//calculate the end of the alignment from b
	end =  bam_calend(&b->core,  bam1_cigar(b));
	end2 = b->core.isize + b->core.pos - 1;
	if(b->core.flag & BAM_FPROPER_PAIR){
		//printf("proper pair qname %s pos1 %d end1 %d pos2 %d end2 %d isize %d \n", bam1_qname(b), b->core.pos, end, b->core.mpos, end2, b->core.isize);
		string name = bam1_qname(b);		
		//check to see if a read has been seen	
		if(buf->find(name)==buf->end()){
			//add alignment
			pair=new pair_t();
			pair->first=NULL;
			pair->second=NULL;
			(*buf)[name]=pair;
			(*plist).push_back(pair);
			
		}
		if(b->core.flag & BAM_FREAD1){
		        (*buf)[name]->first = new bam1_t();	
			(*buf)[name]->first = bam_copy1((*buf)[name]->first, b);
		}else{
			(*buf)[name]->second = new bam1_t();
                        (*buf)[name]->second = bam_copy1((*buf)[name]->second, b);
		}
        }else if(b->core.flag & BAM_FPAIRED && ! (b->core.flag & BAM_FUNMAP) && ! (b->core.flag & BAM_FMUNMAP) && !(b->core.flag &  BAM_FSECONDARY) ){ 
		// not proper pair, but both pairs are mapped and alignment is not a secondary alignment.
	               //printf("proper pair qname %s pos1 %d end1 %d pos2 %d end2 %d isize %d \n", bam1_qname(b), b->core.pos, end, b->core.mpos, end2, b->core.isize);
                string name = bam1_qname(b);
                //check to see if a read has been seen
                if(bad_buf->find(name)==bad_buf->end()){
                        //add alignment
                        pair=new pair_t();
                        pair->first=NULL;
                        pair->second=NULL;
                        (*bad_buf)[name]=pair;
                        (*bad_plist).push_back(pair);

                }
                if(b->core.flag & BAM_FREAD1){
                        (*bad_buf)[name]->first = new bam1_t();
                        (*bad_buf)[name]->first = bam_copy1((*bad_buf)[name]->first, b);
                }else{
                        (*bad_buf)[name]->second = new bam1_t();
                        (*bad_buf)[name]->second = bam_copy1((*bad_buf)[name]->second, b);
                }
	
	}
	return 0;
}
void remove_singlets(list<pair_t *> *plist){
	 //remove singlets
        pair_t *pair;
	for( list<pair_t *>::iterator i=plist->begin(); i != plist->end(); ++i){
	        pair= *i;
                if(!(pair->first && pair->second)){
                	plist->erase(i);
                        i--;
                }
	}
}
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
                        (*mean_insert)+=(double)(stop - start + 1);
                        (*counts)++;
                }
        }
        if((*counts)>0.0){
                (*mean_insert)/=(double)(*counts);
        }

}

void calculate_fragment_coverage(list<pair_t *> *plist, unsigned long *pcov, double *mean_insert, unsigned long * counts, unsigned long *fragment_hist , unsigned int length){

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
			fragment_hist[stop-start+1]++;
                        for(unsigned int count=start; count<=stop;count++){

                                pcov[count]++;
                        }
                        (*mean_insert)+=(double)(stop - start + 1);
                        (*counts)++;
                }
        }
	if((*counts)>0.0){
		(*mean_insert)/=(double)(*counts);
	}

}
void calculate_fragment_hist(list<pair_t *> *plist, unsigned long *fragment_hist, unsigned int length){

        pair_t *pair=NULL;
        //calculate coverage, and mean insert length.
        for(  list<pair_t *>::iterator i=plist->begin(); i != plist->end(); ++i){
                pair = *i;
                if(pair->first && pair->second){
                        unsigned int start, stop;
                        calc_span(&start,&stop,pair);
                        //alignments can go past the end of the read.
                        start = (start <= 0)? 0 : start;
                        stop = (stop >= length)? length-1 : stop;
			fragment_hist[stop-start+1]++;
                }
        }
}

/*
creates a model probability model for testing
*/
model_t::model_t(unsigned long size, unsigned long *insert_counts, unsigned long counts,ostream *logout){
        p = new double[size];
        E = new double[size];
        pI = new double[size];
        EI = new unsigned int[size];
        memset(p,0,sizeof(double)*size);
	memset(E,0,sizeof(double)*size);
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
                        *logout << "model_E\t" << i<< "\t" << E[i] << "\n";
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


void calculate_transcript(list<pair_t *> *plist, list<pair_t *> *bad_plist, unsigned int target_len,ostream *logout){
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
        model_t *model = new model_t(target_len, fragment_hist, counts, logout);
        //calculate the minimum probability of any position in the transcript.
        double min_binom_p=1.0;
        unsigned int min_ink=0;
        double p_total=0.0;
        for(unsigned int ink=0;ink<target_len;++ink){
                double P=0.0;
                if(model->p[ink] < 0.0 || model->p[ink] > 1.0){
                        cerr << "p out-of-bounds " << model->p[ink] << " at " << ink << "\n";
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
                (*logout)<< "pcov_P\t" << ink << "\t" << pcov[ink] <<"\n";
        }
        double res = 0;
        if(N>0 || bad_count==0){
                res = cdf_multinomial_P(target_len,N,model->pI,model->EI);
        }else{
                res = -1;
        }
        cout << "len " << target_len <<" frag_c "<< counts << " tot_Nucs " << N << " ave_frag_len " << mean_insert << " pos " << min_ink <<" minPcdf " << min_binom_p << " model_p " << model->p[min_ink] << " model_expec " << model->E[min_ink] << " m_cdf " << res << " pcov " << pcov[min_ink] << " prop " << plist->size() << " imp " << bad_plist->size() << " p_tot " << p_total <<"\n";
        //cout << "\n";
        delete pcov;
        delete fragment_hist;
        delete model;
}

