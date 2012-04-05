#include <stdio.h>
#include <iostream>
#include <string.h>
#define _cplusplus
#include <list>
#include <cstdio>
#include "fragment_coverage.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <ostream>
#include "model.h"
#include "multinomial.h"


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
void calculate_coverage_counts(list<pair_t *> *plist,unsigned long *cov, unsigned long *mcov, double *mean_insert, unsigned long * counts, unsigned long *fragment_hist , unsigned int length){
/* This method is for calculating the md point scan stat, instead calculating the coverage it counts the positions of the midpoints of the fragment.


*/
        pair_t *pair=NULL;

        //process contig calculating paired insert coverage.
        memset(cov,0,sizeof(long)*length);
        memset(cov,0,sizeof(long)*length);
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
			mcov[midpoint]++;
                        for(unsigned int count=start; count<=stop;count++){

                                cov[count]++;
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




double Fp(int k, double psi){
	//Wrappers for the gsl cumulative poisson.
	//The default for values of k <= 0.0 is 1.0
	//The scan stat book defined them to be 0.0
	  if(k<0 ){
                return 0.0;
        }else {
		return gsl_cdf_poisson_P(k,psi);
	}	
}

double pois(int k, double psi){
	//Wrapper for gsl_poison pdf.
	//This is for convenience of coding and testing different numerical methods.
	if(k<0 ){
		return 0.0;
	}else {

		return gsl_ran_poisson_pdf(k,psi);
	}
}

double r_scan_approximation(long k,double psi, double L){

/*  Calculates the probability that the sum of a contiguous subset from random variates X1,X2, ... XN is smaller than the expected minimum order statistic given the expected number of variates to cover w.
 psi=lambda*w
 L=T/w
*/
        double P=0.0;

        double Q2,Q3,A1,A2,A3,A4,Fp_minus_1,Fp_minus_2,Fp_minus_3,pk;
        Fp_minus_3=Fp(k-3,psi); //evaluate Fb as few times as possible.
        Fp_minus_2=Fp_minus_3+pois ( k-2, psi);
        Fp_minus_1=Fp_minus_2+pois ( k-1, psi);

        pk=pois ( k, psi);
        Q2 = pow(Fp_minus_1,2)-(k-1)*pk*pois ( k-2, psi)
                -(k-1-psi)*pk*Fp_minus_3;
        A1=2*pk*Fp_minus_1*((k-1)*Fp_minus_2-psi*Fp_minus_3);
        A2=.5*pow(pk,2)*((k-1)*(k-2)*Fp_minus_3-
                2*(k-2)*psi*Fp(k-4,psi)+
                psi*psi*Fp(k-5,psi));
        A3=0.0;
        int r=1;

        A4=0.0;
        double Fp_rmin1=0.0,Fp_rmin2=0.0,Fp_rmin3=0.0;
        Fp_rmin1 = pois(0,psi);
        A3+=pois ( 2*k-r, psi)*pow(Fp_rmin1,2);
        for( r=2;r<=k-1;r++){
                Fp_rmin3 = Fp_rmin2;
                Fp_rmin2 += pois(r-2,psi);
                Fp_rmin1 += pois(r-1,psi);
                A3+=pois ( 2*k-r, psi)*pow(Fp_rmin1,2);
                A4+=pois ( 2*k-r, psi)*pois ( r, psi)*((r-1)*Fp_rmin2-psi*Fp_rmin3);

        }

        Q3= pow(Fp_minus_1,3)-A1+A2+A3-A4;
        P=1-Q2*pow(Q3/Q2,L-2);
        return P;
}


double large_counts_r_scan_approximation(long k,double lambda, double T, double w){
/*	This method estimates the r-scan statistic quickly and is accurate for large counts.
	For the purposes of the coverage quality statistic.  While the coverage problem is looking for the minimal statistic, the dual maximual problem is well defined.
	k'=N-k
	w'=T-w

	This method calculates the r_scan_approximation for large counts, with only one call to the iterative convergence numerical algorithm Fp, which is faster than the 
	method that is accurate across all input ranges*/
	//return 1.0 - Fp(k-1,lambda*w)*exp(-(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w));
	double x = -(log(Fp(k-1,lambda*w))-(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w));
//        cout << " x " << x << " P " << 1.0 - Fp(k-1,lambda*w)*exp(-(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w)) << " ratio " << (1.0000001 - 0.16667825 * x) / (1 + 0.33321764 * x) <<"\n";
        double p = x * (1.0000001 - 0.16667825 * x) / (1 + 0.33321764 * x);

	//The following is to pretty up the output changing -0 to just 0.
	if (p==0.0){
		p=0.0;
	}
	return p;

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
        P=large_counts_r_scan_approximation(k,lambda,T,w);
	//P=r_scan_approximation(k,w*lambda, T/w);
	
        //cout << " w " << w << " Ew " << Ew_nonhomo << " EW_homo " << Ew_homo << " ET " << ET << " T " << T << " lambda " << lambda << " psi " << psi <<  " L " << L << " P " << P <<  "\n"
        return P;
}

void fast_md_pt_scan_stat(long N, unsigned long *coverage, unsigned long *fragment_hist, unsigned long counts ,double *p,long *start,long *r, double *exp_win_cov, long * win_cov, ostream *logout){
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
		if(fragment_hist[i] != 0){
			double lambda = ((double)fragment_hist[i])/((double)(obsable_end-obsable_start)+1);
			//(logout << " lambda " << lambda << " obsable_start " << obsable_start << " obsable_end " << obsable_end << " fh " << fragment_hist[i] << " i " << i << "\n";
			for(int j=obsable_start;j<=obsable_end;j++){
				lambdas[j-1]+=lambda;
			}
		}
		//calculate the minimum lenght fragment
		 if(min_f==-1 && fragment_hist[i]!=0){
			//cout << " i " << i << " fh " << fragment_hist[i] << "\n";
                        min_f=i;
			min_w=obsable_start;
			max_w=obsable_end;
                }


	}

	


        for(int i=0;i<N;i++){
                total_std_cov+=coverage[i];
                total_prop_cov+=lambdas[i];
        }
        //standardize expected coverage such that the proportions of the bases in the transcript are the maintained,
        //but the total coverage is the same as for the standarrdized coverage.
	double proportions = ((double)total_std_cov)/total_prop_cov;

//	cout << " proportions " << proportions << " cov " << total_std_cov << " prop " << total_prop_cov << "\n";
        for(int i=0;i<N;i++){
                std_exp_cov[i]=lambdas[i]*proportions;
		(*logout) << "cov\t" << i << "\t" << std_exp_cov[i] << "\t" << coverage[i] << "\n";
        }

        //for each window on the transcript


        double min_p=1.0;
        long min_start=0,min_r=1;
        double temp_p=1.0;
	double min_win_exp_cov=0.0;
	long min_win_cov=0;

        //calculate stat for every pair
        long x2;
        long x1;
        long cum_cov;
        double exp_cov;
        long temp_r;
	long w= max_w - min_w + 1;
	//
        for(x1=min_w;x1<=max_w-1;x1++){
                //each time the cumulants will start at x1;
                cum_cov=coverage[x1-1];
                exp_cov=std_exp_cov[x1-1];

                for(x2=x1+1;x2<=max_w;x2++){
                        cum_cov+=coverage[x2-1];
                        exp_cov+=std_exp_cov[x2-1];
                        temp_r=x2-x1+1;
                        if(cum_cov < exp_cov*.5  ){
				temp_p = fast_min_r_scan_calculate(temp_r, w, cum_cov, total_std_cov, exp_cov);
                        }else{
                                temp_p=1.0;
                        }
                        if(temp_p<min_p || (temp_p==min_p && min_r<temp_r)){
                                min_p=temp_p;
                                min_start=x1;
                                min_r=temp_r;
				min_win_exp_cov=exp_cov;
				 min_win_cov=cum_cov;
/*                                cout << "P\t" << min_p << "\tr\t" << min_r << 
					 " x1 " << x1 << " x2 " << x2 << " cum_cov " << cum_cov
                                         <<" exp_cov " << exp_cov <<  " tot_cov " << total_std_cov <<"\n";*/
			}
                        
                }

        }
	
        (*p)=min_p;
        (*start)=min_start;
        (*r)=min_r;
	(*exp_win_cov)=min_win_exp_cov;
	(*win_cov)=min_win_cov;
        delete std_exp_cov;
}
void minimum_coverage_probability( 
	unsigned long *pcov, 
	unsigned long *fragment_hist, 
	double mean_insert,
	unsigned long counts,
	unsigned int target_len,
	unsigned long *total_nucs,
	unsigned long *position,
	double *model_bp,
	double *model_mcdf,
	double *coverage,
	double *exp_cov,
	ostream *logout){


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
                res = 1.0-cdf_multinomial_P(target_len,N,model->pI,model->EI);
		if(res<0.0){// the approximation can have very small negative numbers, which will be interpretted as 0
			res=0.0;
		}
        }else{
                res = -1;
        }
	(*model_bp)=model->p[min_ink];
	(*model_mcdf)=res;
	(*coverage)= pcov[min_ink];
	(*exp_cov)=model->E[min_ink];
	(*position)=min_ink;
	(*total_nucs)=N;
        delete model;
}

void calculate_transcript_shape_coverage(list<pair_t *> *plist, list<pair_t *> *bad_plist, unsigned int target_len,ostream *logout){
        //process contig calculating paired insert coverage.
        unsigned long *raw_coverage=new unsigned long[target_len];
        memset(raw_coverage,0,sizeof(unsigned long)*target_len);
	unsigned long *midpoint_counts=new unsigned long[target_len];
        memset(midpoint_counts,0,sizeof(unsigned long)*target_len);

        //Clean up hash and list for next set of pairs.
        double mean_fragment=0.0;
        unsigned long counts=0;
        //calculate coverage, and mean insert length.
        unsigned long *fragment_hist=new unsigned long[target_len+1];
	memset(fragment_hist,0,sizeof(unsigned long)*(target_len+1));
        calculate_coverage_counts(plist,raw_coverage,midpoint_counts,&mean_fragment,&counts,fragment_hist,target_len);

        double p=0.0;
        long r=target_len;
        long start=0;
	double win_exp_cov=0;
	long win_cov=0;
	unsigned long total_nucs=0;
	unsigned long position=0;
	double model_bp=0;
	double model_mcdf=0;
	double coverage=0;
	double exp_cov=0;
	
	if(counts>0){

		minimum_coverage_probability( raw_coverage, fragment_hist, mean_fragment, counts, target_len, 
			&total_nucs, &position, &model_bp, &model_mcdf, &coverage, &exp_cov, logout);

	        fast_md_pt_scan_stat(target_len,midpoint_counts,fragment_hist,counts,&p,&start,&r,&win_exp_cov,&win_cov, logout);
	}
	//coverage_scan
        //cout << "len\tfrag_c\tmean_frag\ti\tr\tP(k,w)\n";
        cout << "\tT="<< target_len << "\tF=" << counts << "\tav_f=" << mean_fragment << "\tshape:\tx1=" << start  << "\tr="<<r << "\tp=" << p << "\twin_count=" << win_cov << "\texp_count="<< win_exp_cov << "\tscale:\ttot_nuc="<< total_nucs<< "\tpos="<<position<<"\tcov="<<coverage<<"\texp(cov)="<<exp_cov<<"\tbin-p="<<model_bp<<"\tmcdf="<< model_mcdf<<"\n";
        delete raw_coverage;
        delete fragment_hist;

}



