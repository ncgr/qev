#include <stdio.h>
#include <iostream>
#include <string.h>
#define _cplusplus
#include "fragment_coverage.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <ostream>

using namespace std;
/* this function is designed to test the statistics routines.*/



double r_scan_approximation_debug(long k,double psi, double L){

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
	cout << "pow(Fp_minus_1,2) " << pow(Fp_minus_1,2) << "(k-1)*pk*pois ( k-2, psi)-(k-1-psi)*pk*Fp_minus_3 " << (k-1)*pk*pois ( k-2, psi)
                -(k-1-psi)*pk*Fp_minus_3 <<"\n";

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
	cout << " Q2 " << Q2 << " pow(Q3/Q2,L-2) " << pow(Q3/Q2,L-2) << "\n"; 
        P=1-Q2*pow(Q3/Q2,L-2);
        return P;
}
double large_counts_r_scan_approximation_debug(long k,double lambda, double T, double w){
/*      This method estimates the r-scan statistic quickly and is accurate for large counts.
        For the purposes of the coverage quality statistic.  While the coverage problem is looking for the minimal statistic, the dual maximual problem is well defined.
        k'=N-k
        w'=T-w

        This method calculates the r_scan_approximation for large counts, with only one call to the iterative convergence numerical algorithm Fp, which is faster than the
        method that is accurate across all input ranges*/
        //return 1.0 - Fp(k-1,lambda*w)*exp(-(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w));
	cout << " -(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w) " << -(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w) << " pois(k-1,lambda*w) " << pois(k-1,lambda*w) << "\n";

	//return 1.0 - Fp(k-1,lambda*w)*exp(-(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w));
	//return 1.0 - exp(log(Fp(k-1,lambda*w))-(1.0+(w*lambda)/((double)k))*lambda*pois(k-1,lambda*w));
	//return 1.0 - Fp(k-1,lambda*w)*exp(-(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w));
	//return -(log(Fp(k-1,lambda*w))-(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w));
	double x = -(log(Fp(k-1,lambda*w))-(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w));
	cout << " x " << x << " P " << 1.0 - Fp(k-1,lambda*w)*exp(-(((double)k-w*lambda)/(double)k)*lambda*pois(k-1,lambda*w)) << " ratio " << (1.0000001 - 0.16667825 * x) / (1 + 0.33321764 * x) <<"\n";
	return x * (1.0000001 - 0.16667825 * x) / (1 + 0.33321764 * x);

}

double fast_min_r_scan_calculate_debug(long wp, long T, long kp, long N, double exp_cov){
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
        P=large_counts_r_scan_approximation_debug(k,lambda,T,w);
        //P=r_scan_approximation_debug(k,lambda*w, T/w);

        //cout << " w " << w << " Ew " << Ew_nonhomo << " EW_homo " << Ew_homo << " ET " << ET << " T " << T << " lambda " << lambda << " psi " << psi <<  " L " << L << " P " << P <<  "\n"
        return P;
}



int main(int argc, char *argv[])
{
	double p;
/*	// k=4 psi*L=10 1/L=.1  L=10 psi=1  p=.374
	p = r_scan_approximation(4,1.0, 10.0);
	cout << "k=4 psi*L=10 1/L=.1  L=10 psi=1  p=.374 p*="<<p<<"\n";

	// k=5 psi*L=12 1/L=.25  L=4 psi=3  p=.765
	p = r_scan_approximation(5,3.0, 4.0);
	cout << "k=5 psi*L=12 1/L=.25  L=4 psi=3  p=.765 p*="<<p<<"\n";

	// k= 5 psi*L=8 1/L=1/6  L=6 psi=8/6  p=.896
	p = r_scan_approximation(5,8.0/6.0, 6.0);
	cout << "k= 5 psi*L=8 1/L=1/6  L=6 psi=8/6  p=.896 p*="<<p<<"\n";
	
	// k= 41 psi*L=4800 1/L=1/240  L=240 psi=20  p=.066
        p = r_scan_approximation(41,20, 240);
        cout << "k= 41 psi*L=4800 1/L=1/240  L=240 psi=20  p=.066 p*="<<p<<"\n";

	// k= 11 psi*L=33 1/L=1/6  L=6 psi=33/6  p=.381
        p = r_scan_approximation(11,33.0/6.0, 6.0);
        cout << "k= 11 psi*L=33 1/L=1/6  L=6 psi=33/6  p=.381 p*="<<p<<"\n";
	
	// k= 6 psi*L=342 1/L=1/674.51  L=674.51 psi=.5070347  p=.051
        p = r_scan_approximation(6,342.0/674.51, 674.51);
        cout << "k= 6 psi*L=342 1/L=1/674.51  L=674.51 psi=.5070347  p=.051 p*="<<p<<"\n";

*/
	p = fast_min_r_scan_calculate_debug(1078, 12174-308, 155, 3498, 317.754);
	cout << "test p " << p << " >BPA_1  len     12174   counts  3498    mean_frag       308.244 start   10565   r       1078    p       0.00214578      win_cov 155     win_exp_cov     317.754\n";

	p = fast_min_r_scan_calculate_debug(2377, 6525-307, 604, 3300, 1208.25);
        cout << "test p " << p << " >BPA_26         len     6525    counts  3300    mean_frag       307.24  start   4096    r       2377    p       0       win_cov 604     win_exp_cov     1208.25\n";



}

/*
>BPA_1  len     12174   counts  3498    mean_frag       308.244 start   10565   r       1078    p       0.00214578      win_cov 155     win_exp_cov     317.754
>BPA_2  len     11058   counts  3039    mean_frag       304.425 start   329     r       1225    p       0.000491792     win_cov 173     win_exp_cov     346.145
>BPA_3  len     11424   counts  4182    mean_frag       305.432 start   167     r       881     p       0.00383248      win_cov 164     win_exp_cov     331.001
>BPA_4  len     8153    counts  6471    mean_frag       303.583 start   3468    r       286     p       0.0567519       win_cov 110     win_exp_cov     235.736
>BPA_5  len     10711   counts  2494    mean_frag       305.455 start   8775    r       622     p       0.0623968       win_cov 74      win_exp_cov     149.061
>BPA_6  len     7807    counts  7255    mean_frag       306.822 start   3884    r       824     p       9.8573e-08      win_cov 374     win_exp_cov     796.922
>BPA_7  len     13365   counts  3518    mean_frag       303.697 start   10816   r       2325    p       4.40115e-09     win_cov 311     win_exp_cov     626.161
>BPA_8  len     7434    counts  3456    mean_frag       305.158 start   152     r       874     p       7.61709e-05     win_cov 209     win_exp_cov     420.402
>BPA_9  len     7552    counts  4239    mean_frag       302.803 start   5649    r       388     p       0.00358975      win_cov 55      win_exp_cov     226.844
>BPA_10         len     8027    counts  4241    mean_frag       307.058 start   7012    r       627     p       0.00210319      win_cov 164     win_exp_cov     344.387
>BPA_11         len     7784    counts  4179    mean_frag       306.58  start   4350    r       830     p       8.01645e-05     win_cov 231     win_exp_cov     463.792
>BPA_12         len     12144   counts  2736    mean_frag       303.779 start   10581   r       1429    p       0.000396462     win_cov 159     win_exp_cov     325.93
>BPA_13         len     7134    counts  3651    mean_frag       305.145 start   3318    r       454     p       0.0194085       win_cov 121     win_exp_cov     242.681
>BPA_14         len     6542    counts  10977   mean_frag       303.95  start   1257    r       237     p       0.0150722       win_cov 193     win_exp_cov     416.957
>BPA_15         len     7482    counts  2940    mean_frag       307.421 start   1109    r       382     p       0.0181685       win_cov 45      win_exp_cov     156.507
>BPA_16         len     7987    counts  3239    mean_frag       305.909 start   6245    r       1201    p       3.04191e-08     win_cov 218     win_exp_cov     506.356
>BPA_17         len     7829    counts  2281    mean_frag       306.676 start   2184    r       377     p       0.103667        win_cov 55      win_exp_cov     114.298
>BPA_18         len     8035    counts  2138    mean_frag       305.327 start   4341    r       332     p       0.13799 win_cov 42      win_exp_cov     91.8141
>BPA_19         len     12091   counts  5457    mean_frag       304.362 start   405     r       1090    p       0.000179557     win_cov 251     win_exp_cov     504.591
>BPA_20         len     7093    counts  3103    mean_frag       306.882 start   5474    r       682     p       7.08957e-05     win_cov 108     win_exp_cov     311.788
>BPA_21         len     7253    counts  5986    mean_frag       305.485 start   199     r       846     p       2.92489e-07     win_cov 362     win_exp_cov     728.776
>BPA_22         len     8673    counts  2217    mean_frag       309.19  start   259     r       599     p       0.0409937       win_cov 79      win_exp_cov     158.751
>BPA_23         len     6098    counts  2963    mean_frag       306.444 start   1414    r       343     p       0.02974 win_cov 75      win_exp_cov     175.443
>BPA_24         len     6234    counts  4652    mean_frag       304.921 start   51      r       882     p       1.84488e-09     win_cov 231     win_exp_cov     612.006
>BPA_25         len     6391    counts  5390    mean_frag       303.852 start   3147    r       584     p       0.000113487     win_cov 257     win_exp_cov     517.005
>BPA_26         len     6525    counts  3300    mean_frag       307.24  start   4096    r       2377    p       0       win_cov 604     win_exp_cov     1208.25


*/
