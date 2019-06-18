#include "model/model.h"
#include <stdio.h>

int main() {
    // misc
    int i=0;
    int j=0;
    int k=0;
    // simulation parameters
    int dataPoints=(int)1e1;
    double dt=1e-7;
    // model parameters
    int nComps=109;
    int nTypes=3;
    int nElems=nTypes*nComps;
    int* epsilon=(int *)malloc(nTypes*sizeof(int));
    epsilon[0]= 250;// equivalent to 2.5
    epsilon[1]=   1;// equivalent to 0.01
    epsilon[2]=5000;// equivalent to 50
    int capacity=600;
    int* curState=(int *)malloc(nElems*sizeof(int));// initial condition is uniform distribution
    for(i=0;i<nElems;i+=1) {
        if(i<nComps) {
            curState[i]=279;
        } else if(i<2*nComps) {
            curState[i]=81;
        } else {
            curState[i]=139;
        }
    }
    // rng parameters
    int rng_seed=123;
    // output
    int* output=(int *)malloc(dataPoints*nElems*sizeof(int));
    
    // run the model
    getSeries(dataPoints, dt, nComps, nTypes, epsilon, capacity, curState, rng_seed, output);

    // output the results
    for(i=0;i<dataPoints;i+=1) {
        for(j=0;j<nElems;j+=1) {
            k=i*nElems+j;
            printf("%d ",output[k]);
        }
        printf("\n");
    }
}
