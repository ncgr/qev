#include <math.h>
#include "nrutil.h"
#include "num_routines.h"

float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}




float erfcc(float x)
{
	float t,z,ans;

	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	return x >= 0.0 ? ans : 2.0-ans;
}
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort2(unsigned long n, float arr[], float brr[])
{
	unsigned long i,ir=n,j,k,l=1,*istack;
	int jstack=0;
	float a,b,temp;

	istack=lvector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				b=brr[j];
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					brr[i+1]=brr[i];
				}
				arr[i+1]=a;
				brr[i+1]=b;
			}
			if (!jstack) {
				free_lvector(istack,1,NSTACK);
				return;
			}
			ir=istack[jstack];
			l=istack[jstack-1];
			jstack -= 2;
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			SWAP(brr[k],brr[l+1])
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
				SWAP(brr[l],brr[ir])
			}
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
				SWAP(brr[l+1],brr[ir])
			}
			if (arr[l] > arr[l+1]) {
				SWAP(arr[l],arr[l+1])
				SWAP(brr[l],brr[l+1])
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			b=brr[l+1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
				SWAP(brr[i],brr[j])
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			brr[l+1]=brr[j];
			brr[j]=b;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort2.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}

float betai(float a, float b, float x)
{
	float betacf(float a, float b, float x);
	float gammln(float xx);
	void nrerror(char error_text[]);
	float bt;

	if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

float betacf(float a, float b, float x)
{
	void nrerror(char error_text[]);
	int m,m2;
	float aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
	return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN

void spear(double data1[], double data2[],  unsigned long n, double *d, double *zd,
double *probd, double *rs, double *probrs)
//Given two data arrays, data1[1..n] and data2[1..n], this routine returns their sum-squared
//dirence of ranks as D, the number of standard deviations by which D deviates from its nullhypothesis
//expected value as zd, the two-sided signi^Lcance level of this deviation as probd,
//Spearman's rank correlation rs as rs, and the two-sided signi^Lcance level of its deviation from
//zero as probrs. The external routines crank (below) and sort2 (x8.2) are used. A small value
//of either probd or probrs indicates a signi^Lcant correlation (rs positive) or anticorrelation
//(rs negative).
{
        float betai(float a, float b, float x);
        void crank(unsigned long n, float w[], float *s);
        float erfcc(float x);
        void sort2(unsigned long n, float arr[], float brr[]);
        unsigned long j;
        float vard,t,sg,sf,fac,en3n,en,df,aved,*wksp1,*wksp2;
        wksp1=vector(1,n);
        wksp2=vector(1,n);
        for (j=1;j<=n;j++) {
                wksp1[j]=data1[j];
                wksp2[j]=data2[j];
        }
        sort2(n,wksp1,wksp2); //Sort each of the data arrays, and convert the entries to
                                //ranks. The values sf and sg return the sums
                        //P(f3k.fk) and P(g3m. gm), respectively.
        crank(n,wksp1,&sf);
        sort2(n,wksp2,wksp1);
        crank(n,wksp2,&sg);
        *d=0.0;
        for (j=1;j<=n;j++) //Sum the squared dirence of ranks.
                *d += SQR(wksp1[j]-wksp2[j]);
        
        en=n;
        en3n=en*en*en-en;
        aved=en3n/6.0-(sf+sg)/12.0;// Expectation value of D,
        fac=(1.0-sf/en3n)*(1.0-sg/en3n);
        vard=((en-1.0)*en*en*SQR(en+1.0)/36.0)*fac; //and variance of D give
        *zd=(*d-aved)/sqrt(vard); //number of standard devia-
        *probd=erfcc(fabs(*zd)/1.4142136);// tions and signi^Lcance.
        *rs=(1.0-(6.0/en3n)*(*d+(sf+sg)/12.0))/sqrt(fac);// Rank correlation coecient,
        fac=(*rs+1.0)*(1.0-(*rs));
        if (fac > 0.0) {
                t=(*rs)*sqrt((en-2.0)/fac);// and its t value,
                df=en-2.0;
                *probrs=betai(0.5*df,0.5,df/(df+t*t));// give its signi^Lcance.
        } else
                *probrs=0.0;
        free_vector(wksp2,1,n);
        free_vector(wksp1,1,n);
        
}
void crank(unsigned long n, float w[], float *s)
//Given a sorted array w[1..n], replaces the elements by their rank, including midranking of ties,
//and returns as s the sum of f3 . f, where f is the number of elements in each tie.
{
        unsigned long j=1,ji,jt;
        float t,rank;
        *s=0.0;
        while (j < n) {
                if (w[j+1] != w[j]) {// Not a tie.
                        w[j]=j;
                        ++j;
                } else {// A tie:
                        for (jt=j+1;jt<=n && w[jt]==w[j];jt++);// How far does it go?
                        rank=0.5*(j+jt-1); //This is the mean rank of the tie,
                        for (ji=j;ji<=(jt-1);ji++) w[ji]=rank;// so enter it into all the tied
                        t=jt-j;// entries,
                        *s += t*t*t-t;// and update s.
                        j=jt;
                }
        }
        if (j == n) w[n]=n;// If the last element was not tied, this is its rank.
}

#include <math.h>
#define TINY 1.0e-20 //Will regularize the unusual case of complete correlation.
void pearsn(double x[], double y[], unsigned long n, double *r, double *prob, double *z)
/*Given two arrays x[1..n] and y[1..n], this routine computes their correlation coecientcient
r (returned as r), the signicance level at which the null hypothesis of zero correlation is
disproved (prob whose small value indicates a signicant correlation), and Fisher's z (returned
as z), whose value can be used in further statistical tests as described above.
*/
{
	float betai(float a, float b, float x);
	float erfcc(float x);
	unsigned long j;
	double yt,xt,t,df;
	double syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;
	for (j=1;j<=n;j++) { //Find the means.
		ax += x[j];
		ay += y[j];
	}
	ax /= n;
	ay /= n;
	for (j=1;j<=n;j++) { //Compute the correlation coecientcient.
		xt=x[j]-ax;
		yt=y[j]-ay;
		sxx += xt*xt;
		syy += yt*yt;
		sxy += xt*yt;
	}
	*r=sxy/(sqrt(sxx*syy)+TINY);
	*z=0.5*log((1.0+(*r)+TINY)/(1.0-(*r)+TINY));//Fisher's z transformation.
	df=n-2;
	t=(*r)*sqrt(df/((1.0-(*r)+TINY)*(1.0+(*r)+TINY))); //Equation (14.5.5).
	*prob=betai(0.5*df,0.5,df/(df+t*t)); //Student's t probability.
	/* *prob=erfcc(fabs((*z)*sqrt(n-1.0))/1.4142136) */
	//For large n, this easier computation of prob, 
	//using the short routine erfcc, would give approximately
	//the same value.
}




#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7

void gser(float *gamser, float a, float x, float *gln)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	int n;
	float sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(float *gammcf, float a, float x, float *gln)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	int i;
	float an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */
float gammq(float a, float x)
{
	void gcf(float *gammcf, float a, float x, float *gln);
	void gser(float *gamser, float a, float x, float *gln);
	void nrerror(char error_text[]);
	float gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */



