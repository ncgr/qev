#include <stdio.h>
#define _cplusplus
#include<ostream>
using namespace std;



class model_t{
        public:
        double *p;
        double *E;
        double *pI;
        unsigned int *EI;
        model_t(unsigned long size, unsigned long *insert_counts, unsigned long counts,ostream *logout);
        ~model_t();
};


