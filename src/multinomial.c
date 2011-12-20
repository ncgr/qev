#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

/* Computes the (log) CDF of a multinomial distribution with parameters p and N
 * i.e. P(X1<=n1,X2<=n2,...,XK<=nK)
 * Uses a approximation from 
 * "A Representation for Multinomial Cumulative Distribution Functions", 
 * Bruce Levin, The Annals of Statistics, v.9, n.5, pp.1123-1126, 1981
 */

double
cdf_multinomial_lnP (const size_t K, const unsigned int N,
                           const double p[], const unsigned int n[])
{
  size_t k;
  double log_cdf = 0.0;
  double s;
  double gamma1=0.0, gamma2=0.0, sum_s2=0.0, sum_mu=0.0;
  double mu, s2, sp, pcdf, mr, mf2, mf3, mf4, mu2, mu3, mu4, x, x2, PWN;
  double phi=0.0;
  s = (double) N;
  log_cdf = -log(gsl_ran_poisson_pdf(N,s));

  /* This is the P(W=N) bit */
  for(k=0; k<K; k++) {
    sp = s*p[k];

    //    pcdf = gsl_cdf_poisson_P(n[k],sp); original
    pcdf = 1-gsl_cdf_poisson_P(n[k],sp)+gsl_ran_poisson_pdf(n[k],sp); // rsk:modified to calculate other tail
    log_cdf += log(pcdf);

   //mu = sp*(1-gsl_ran_poisson_pdf(n[k],sp)/pcdf); original
    if(n[k]<=1.0){
	mu=sp;
	s2=mu;
    }else{
    	mu = sp*(1+(gsl_ran_poisson_pdf(n[k]-1,sp))/pcdf); //rsk:modified to change truncated mean to reflect <= to => change
    //s2 = mu-(n[k]-mu)*(sp-mu)  modified to refelct changes in mu
    
    	phi=sp*sp*(1+(gsl_ran_poisson_pdf(n[k]-2,sp)+gsl_ran_poisson_pdf(n[k]-1,sp))/pcdf);
    	s2 = phi+mu-mu*mu;
    }
    /* factorial moments? */
    mr = n[k];
    mf2= sp*mu-mr*(sp-mu);

    mr *= n[k]-1;
    mf3 = sp*mf2-mr*(sp-mu);
    
    mr *= n[k]-2;
    mf4 = sp*mf3-mr*(sp-mu);

    /* Central Moments */
    mu2 = mf2+mu*(1-mu);
    mu3 = mf3+mf2*(3-3*mu)+mu*(1+mu*(-3+2*mu));
    mu4 = mf4+mf3*(6-4*mu)+mf2*(7+mu*(-12+6*mu))+mu*(1+mu*(-4+mu*(6-3*mu)));

    /* accumulate coef skewness and excess */
    gamma1 += mu3;
    gamma2 += mu4-3*s2*s2;
    sum_mu += mu;
    sum_s2 += s2;
  }
  
  sp = sqrt(sum_s2);
  gamma1 /= sum_s2*sp;
  gamma2 /= sum_s2*sum_s2;

  x = (N-sum_mu)/sp;
  x2 = x*x;
  PWN = -x2/2
    +log(1+gamma1/6*x*(x2-3)+gamma2/24*(x2*x2-6*x2+3)+
	 gamma1*gamma1/72*(((x2-15)*x2+45)*x2-15))
    -log(2*M_PI)/2 -log(sp);
  //printf("sum_mu %g sum_s2 %g log_cdf %g x %g f1 %g f2 %g f3 %g\n",sum_mu,sum_s2,log_cdf,x,1+gamma1/6*x*(x2-3),gamma2/24*(x2*x2-6*x2+3), gamma1*gamma1/72*(((x2-15)*x2+45)*x2-15));
  log_cdf += PWN;

  return log_cdf;
}

double
cdf_multinomial_P (const size_t K, const unsigned int N,
                         const double p[], const unsigned int n[])
{
  return exp (cdf_multinomial_lnP (K, N, p, n));
}

/* Computes the Bonferroni-Mallows bounds for the CDF */

void
cdf_multinomial_BMboundsP (const size_t K, const unsigned int N,
			   const double p[], const unsigned int n[], double *lb, double *ub)
{
  double LB = 1.0, UB = 1.0;
  size_t k;
  
  for(k=0; k<K; k++) {
    LB -= gsl_cdf_binomial_Q(n[k],p[k],N);
    UB *= gsl_cdf_binomial_P(n[k],p[k],N); 
  }

  *lb = LB;
  *ub = UB;
}
  

  
/* Special case with all p equal and all n equal
 * ei. P(max_k X_k <= n) when p_k = p for all k
 */

double
cdf_multinomial_lnPmax (const size_t K, const unsigned int N,
			const unsigned int n)
{
  double log_cdf = 0.0;
  double p, s;
  double gamma1, gamma2;
  double mu, s2, sp, pcdf, mr, mf2, mf3, mf4, mu2, mu3, mu4, x, x2, PWN;
  
  p = 1/((double) K);
  s = (double) N;

  log_cdf = -log(gsl_ran_poisson_pdf(N,s));

  sp = s*p;
  pcdf = gsl_cdf_poisson_P(n,sp);
  log_cdf += K*log(pcdf);

  /* This is the P(W=N) bit */
  mu = sp*(1-gsl_ran_poisson_pdf(n,sp)/pcdf);
  s2 = mu-(n-mu)*(sp-mu);

  /* factorial moments? */
  mr = n;
  mf2= sp*mu-mr*(sp-mu);

  mr *= n-1;
  mf3 = sp*mf2-mr*(sp-mu);
    
  mr *= n-2;
  mf4 = sp*mf3-mr*(sp-mu);

    /* Central Moments */
  mu2 = mf2+mu*(1-mu);
  mu3 = mf3+mf2*(3-3*mu)+mu*(1+mu*(-3+2*mu));
  mu4 = mf4+mf3*(6-4*mu)+mf2*(7+mu*(-12+6*mu))+mu*(1+mu*(-4+mu*(6-3*mu)));

    /* accumulate coef skewness and excess */
  sp = sqrt(K*s2);
  gamma1 = mu3/(sp*s2);
  gamma2 = (mu4-3*s2*s2)/(s2*s2*K);

  x = (N-K*mu)/sp;
  x2 = x*x;
  PWN = -x2/2
    +log(1+gamma1/6*x*(x2-3)+gamma2/24*(x2*x2-6*x2+3)+
	 gamma1*gamma1/72*(((x2-15)*x2+45)*x2-15))
    -log(2*M_PI)/2 -log(sp);

  log_cdf += PWN;

  return log_cdf;
}

double
cdf_multinomial_Pmax (const size_t K, const unsigned int N,
		      const unsigned int n)
{
  return exp (cdf_multinomial_lnPmax (K, N, n));
}
