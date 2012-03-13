#include <stdio.h>
#define _cplusplus
#include <list>
#include<ostream>
#include "bam_coverage.h"


using namespace std;


//call back for fetching pairs

double r_scan_approximation(long k,double psi, double L);
void calculate_transcript_shape_coverage(list<pair_t *> *plist, list<pair_t *> *bad_plist, unsigned int target_len,ostream *logout);

double Fp(int k, double psi);

double pois(int k, double psi);


