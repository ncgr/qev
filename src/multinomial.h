extern double
cdf_multinomial_lnP (const size_t K, const unsigned int N,
		     const double p[], const unsigned int n[]);

extern double
cdf_multinomial_P (const size_t K, const unsigned int N,
		   const double p[], const unsigned int n[]);

extern void
cdf_multinomial_BMboundsP (const size_t K, const unsigned int N,
			   const double p[], const unsigned int n[], double *lb, double *ub);


extern double
cdf_multinomial_lnPmax (const size_t K, const unsigned int N,
			const unsigned int n);

extern double
cdf_multinomial_Pmax (const size_t K, const unsigned int N,
		      const unsigned int n);
