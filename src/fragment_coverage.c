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
#include <gsl/gsl_randist.h>
#include <math.h>
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
void calculate_mid_pt_coverage(list<pair_t *> *plist, unsigned long *pcov, double *mean_insert, unsigned long * counts, unsigned long *fragment_hist , unsigned int length){
/* This method is for calculating the md point scan stat, instead calculating the coverage it counts the positions of the midpoints of the fragment.


*/
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

			//start and stop are zero based
                        start = (start <= 0)? 0 : start;
                        stop = (stop >= length)? length-1 : stop;
                        fragment_hist[stop-start+1]++;
			//count only the mid-point of each fragment
			int midpoint= (long)(.5 + ((double)start+(double)stop)/2.0);
			pcov[midpoint]++;
                        /*for(unsigned int count=start; count<=stop;count++){

                                pcov[count]++;
                        }*/
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




double Fp(int k, double psi){
	//Wrappers for the gsl cumulative poisson.
	//The default for values of k <= 0.0 is 1.0
	//The scan stat book defined them to be 0.0
	return (k<0)?0:gsl_cdf_poisson_P(k,psi);
}

double pois(int k, double psi){
	//Wrapper for gsl_poison pdf.
	//This is for convenience of coding and testing different numerical methods.
	return (k<0)?0:gsl_ran_poisson_pdf(k,psi);
}

double large_counts_r_scan_approximation(long k,double lambda, double T, double w){
/*	This method estimates the r-scan statistic quickly and is accurate for large counts.
	For the purposes of the coverage quality statistic.  While the coverage problem is looking for the minimal statistic, the dual maximual problem is well defined.
	k'=N-k
	w'=T-w

	This method calculates the r_scan_approximation for large counts, with only one call to the iterative convergence numerical algorithm Fp, which is faster than the 
	method that is accurate across all input ranges*/
	return 1.0 - Fp(k-1,lambda*w)*exp(-(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w));
}

double fast_min_r_scan_calculate(long wp, long T, long kp, long N, double exp_cov){
        /*
                Glaz Naus Walenstein, Scan Statistics, pg 325
                "There is a relationship between the maximum and minimum rth order scan"
                Such that the probability of observing a minimum coverage on an interval is equivalent to observing a maximum of coverage in the complement.
        */
        //calculate the complement window w and counts k.
        double w=T-wp;
        long k=N-kp;

        //calculate parameters for scan statistic
        double lambda=((double)N-exp_cov)/(double)w;  //The expected rate on the window.
        //double psi=w*lambda;
        //double L=((double)T)/(double)w;
        //calculate extremal scan statistic.
        double P=1.0;
        
	//if kp is not small say half of all counts then it can't possibly be significant.
	if(lambda>0){
	        P=large_counts_r_scan_approximation(k,lambda,T,w);
	}
        //cout << " w " << w << " Ew " << Ew_nonhomo << " EW_homo " << Ew_homo << " ET " << ET << " T " << T << " lambda " << lambda << " psi " << psi <<  " L " << L << " P " << P <<  "\n";
        return P;
}

void fast_md_pt_scan_stat(long N, unsigned long *coverage, unsigned long *fragment_hist, unsigned long counts ,double *p,long *start,long *r){
        /*finds the most significant window and its significance.
        p,start, and r are return values.

        p is the significance of the minimum coverage window
        start is the start position of the minimum coverage window
        r is the size of the minimum coverage window.

        N transcript length
        coverage is an array of observed coverages of length N
        mean_fragment_length is the mean of the fragments observed in coverage.
        fragments is the number of fragments used in coverage estimate.
        total_coverage is the total observed coverage.
        */

        double *std_exp_cov = new double[N];
        double total_prop_cov=0.0;
        long total_std_cov=0;
	long min_w=-1,max_w=-1,min_f=-1;

	//calculate the mid point lambdas from the fragment distribution
	double *lambdas=new double[N];
	memset(lambdas,0,sizeof(double)*N);
	
	for(int i=1; i<=N;i++){
		//for each fragment length calculate the rates for that fragment, and add to the relevant positions.
		long obsable_start,obsable_end;
		if(i%2==0){
			obsable_start=i/2;
		}else{
			obsable_start=(long)floor(((double)i)/2.0)+1;

		}
		obsable_end=N-obsable_start+1;
		for(int j=obsable_start;j<obsable_end;j++){
			lambdas[j-1]+=((double)fragment_hist[i])/((double)(obsable_end-obsable_start)+1);
		}

		//calculate the minimum lenght fragment
		 if(min_f==-1 && fragment_hist[i]!=0){
			cout << " i " << i << " fh " << fragment_hist[i] << "\n";
                        min_f=i;
			min_w=obsable_start;
			max_w=obsable_end;
                }


	}

	


        for(int i=0;i<=N;i++){
                total_std_cov+=coverage[i];
                total_prop_cov+=lambdas[i];
        }
        //standardize expected coverage such that the proportions of the bases in the transcript are the maintained,
        //but the total coverage is the same as for the standarrdized coverage.
        for(int i=0;i<N;i++){
                std_exp_cov[i]=lambdas[i]*total_std_cov/total_prop_cov;
        }

        //for each window on the transcript


        double min_p=1.0;
        long min_start=0,min_r=1;
        double temp_p=1.0;

        //calculate stat for every pair
        long x2;
        long x1;
        long cum_cov;
        double exp_cov;
        long temp_r;
        for(x1=min_w;x1<=max_w-1;x1++){
                //each time the cumulants will start at x1;
                cum_cov=coverage[x1-1];
                exp_cov=std_exp_cov[x1-1];

                for(x2=x1+1;x2<=max_w /*&& x2-x1+1<=1000*/;x2++){
                        cum_cov+=coverage[x2-1];
                        exp_cov+=std_exp_cov[x2-1];
                        temp_r=x2-x1+1;
                        if(cum_cov < .5 * exp_cov){

				temp_p = fast_min_r_scan_calculate(temp_r, max_w - min_w + 1, cum_cov, total_std_cov, exp_cov);
                        }else{
                                temp_p=1.0;
                        }
                        if(temp_p<min_p || (temp_p==min_p && min_r<temp_r)){
                                min_p=temp_p;
                                min_start=x1;
                                min_r=temp_r;
			
                                cout << "P\t" << min_p << " tp " << temp_p<< "\tr\t" << min_r << "\tmin_start\t" << min_start << 
					 " x1 " << x1 << " x2 " << x2 << " cum_std_cov " << cum_cov
                                         << " total_std_cov " << total_std_cov <<" exp_cov " << exp_cov << "\n";
			}
                        
                }

        }
        (*p)=min_p;
        (*start)=min_start;
        (*r)=min_r;
        delete std_exp_cov;
}


void calculate_transcript_scan_stat_mid_pt(list<pair_t *> *plist, list<pair_t *> *bad_plist, unsigned int target_len,ostream *logout){
        //process contig calculating paired insert coverage.
        unsigned long *raw_coverage=new unsigned long[target_len];
        memset(raw_coverage,0,sizeof(unsigned long)*target_len);
        //Clean up hash and list for next set of pairs.
        double mean_fragment=0.0;
        unsigned long counts=0;
        //calculate coverage, and mean insert length.
        unsigned long *fragment_hist=new unsigned long[target_len+1];
	memset(fragment_hist,0,sizeof(unsigned long)*target_len+1);
        calculate_mid_pt_coverage(plist,raw_coverage,&mean_fragment,&counts,fragment_hist,target_len);

        double p=0.0;
        long r=0;
        long start=0;

        fast_md_pt_scan_stat(target_len,raw_coverage,fragment_hist,counts,&p,&start,&r);
        cout << "len\tfrag_c\tmean_frag\ti\tr\tP(k,w)\n";
        cout << target_len << "\t" << counts << "\t" << mean_fragment << "\t" << start  << "\t"<<r << "\t" << p << "\n";
        delete raw_coverage;
        delete fragment_hist;


}


