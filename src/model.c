#define _cplusplus
#include <math.h>
#include <ostream>
#include "model.h"

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


