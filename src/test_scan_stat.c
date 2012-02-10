#include <stdio.h>
#include <iostream>
#include <string.h>
#define _cplusplus
#include "fragment_coverage.h"

using namespace std;
/* this function is designed to test the statistics routines.*/

int main(int argc, char *argv[])
{
	double p;
	// k=4 psi*L=10 1/L=.1  L=10 psi=1  p=.374
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


}

