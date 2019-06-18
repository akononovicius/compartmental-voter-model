#include "model.h"

/*
 * Transition rate formula incorporating all restrictions.
 *
 * Input:
 *      nComps   - number of compartments
 *      i        - state id from which the transition would occur. Note that
 *                 state id includes information about both compartment and
 *                 agent type:
 *                     sid = type*nComps + compartment
 *      j        - state id to which the transition would occur. See info about
 *                 i.
 *      epsilon  - rate of individual transition (depends on type; see
 *                 documentation of getSeries for some caveats)
 *      capacity - capacity of compartments
 *      curState - current occupation of the states. Array indexed by state id
 *                 (see info about i).
 *      popState - current occupation of the compartments
 */
long transitionRate(int nComps, int i, int j, int* epsilon, int capacity,
                   int* curState, int* popState) {
    // self looping transition are forbidden
    if(i==j) {
        return 0;
    }
    // it is forbidden to change type
    int typeFrom=(int)floor(i/nComps);
    int typeTo=(int)floor(j/nComps);
    if(typeFrom!=typeTo) {
        return 0;
    }
    // it is forbidden to go over capacity
    int popTo=popState[j % nComps];
    if(popTo>=capacity) {
        return 0;
    }
    // calculate rate based on homophily formula:
    //      la[ i->j | k ] = x[i]*(e[k] + x[j]);
    // x[j] was multiplicated by 100, because epsilons are already multiplied by
    // 100 (this is needed to keep the balance between individual transitions
    // and homophily induced transitions; model simply speeds up by a factor of
    // 100)
    return curState[i]*(epsilon[typeFrom]+100*curState[j]);
}

/*
 * Initializes transition rate matrix, its row sums and total sum.
 *
 * Input:
 *      nComps   - number of compartments
 *      nTypes   - number of types
 *      epsilon  - rate of individual transition (depends on type; see
 *                 documentation of getSeries for some caveats)
 *      capacity - capacity of compartments
 *      curState - current occupation of the states. Array indexed by state id.
 *                 Note that state id includes information about both
 *                 compartment and agent type:
 *                     sid = type*nComps + compartment
 *      popState - current occupation of the compartments
 * Output:
 *      matrix   - transition rate matrix
 *      rowSum   - sums of each row of the transition rate matrix
 *      totalSum - total sum of all elements in the matrix
 */
void getTransitionRateMatrix(int nComps, int nTypes, int* epsilon, int capacity,
                             int* curState, int* popState, long* matrix,
                             long* rowSum, long* totalSum) {
    int i=0;
    int j=0;
    int k=0;
    int nElems=nComps*nTypes;

    // fill the transition rate matrix
    //  and calculate the total and row sums (which will be used for faster
    //  searching within the matrix)
    *totalSum=0;
    for(i=0;i<nElems;i+=1) {
        rowSum[i]=0;
        for(j=0;j<nElems;j+=1) {
            k=i*nElems+j;// position in the flattened matrix
            matrix[k]=transitionRate(nComps,i,j,epsilon,capacity,curState,
                                     popState);
            rowSum[i]+=matrix[k];
        }
        *totalSum=(*totalSum)+rowSum[i];
    }
}

/*
 * Updates transition rate matrix, its row sums and total sum.
 * This function is useful as it is faster to update the matrix, instead of
 * recalculating it from scratch.
 *
 * Input:
 *      compId   - compartment id to update
 *      nComps   - number of compartments
 *      nTypes   - number of types
 *      epsilon  - rate of individual transition (depends on type; see
 *                 documentation of getSeries for some caveats)
 *      capacity - capacity of compartments
 *      curState - current occupation of the states. Array indexed by state id.
 *                 Note that state id includes information about both
 *                 compartment and agent type:
 *                     sid = type*nComps + compartment
 *      popState - current occupation of the compartments
 * Output:
 *      matrix   - transition rate matrix
 *      rowSum   - sums of each row of the transition rate matrix
 *      totalSum - total sum of all elements in the matrix
 */
void updateTransitionRateMatrix(int compId, int nComps, int nTypes,
                                int* epsilon, int capacity, int* curState,
                                int* popState, long* matrix, long* rowSum,
                                long* totalSum) {
    int i=0;
    int j=0;
    int k=0;
    int m=0;
    int nElems=nComps*nTypes;

    // transition rates should be updated for all agents types within the
    //  compartment
    for(i=0;i<nTypes;i+=1) {
        // transition rates should be update for a whole row and column
        //  corresponding to the compartment and agent type
        for(j=0;j<nElems;j+=1) {
            m=compId+i*nComps;

            // column update (transitions to compartment, which is updated)
            k=j*nElems+m;
            rowSum[j]-=matrix[k];
            *totalSum-=matrix[k];
            matrix[k]=transitionRate(nComps,j,m,epsilon,capacity,curState,
                                     popState);
            rowSum[j]+=matrix[k];
            *totalSum+=matrix[k];

            // row update (transitions from compartment, which is updated)
            k=m*nElems+j;
            rowSum[m]-=matrix[k];
            *totalSum-=matrix[k];
            matrix[k]=transitionRate(nComps,m,j,epsilon,capacity,curState,
                                     popState);
            rowSum[m]+=matrix[k];
            *totalSum+=matrix[k];
        }
    }
}

/*
 * Updates the current model state by doing a single step using Gillespie
 * method.
 *
 * Input:
 *      nComps   - number of compartments
 *      nTypes   - number of types
 *      epsilon  - rate of individual transition (depends on type; see
 *                 documentation of getSeries for some caveats)
 *      capacity - capacity of compartments
 *      curState - current occupation of the states. Array indexed by state id.
 *                 Note that state id includes information about both
 *                 compartment and agent type:
 *                     sid = type*nComps + compartment
                   (also output)
 *      popState - current occupation of the compartments (also output)
 *      clock    - number of time units passed (also output)
 *      matrix   - array allocated for the transition matrix
 *      rowSum   - array allocated for the row sums of the transition matrix
 *      totalSum - total sum of the transition matrix
 *      rexp     - random number drawn from exponential distribution
 *      runi     - random number drawn from uniform distribution
 * Output:
 *      curState - see input
 *      popState - see input
 *      clock    - see input
 */
void gillespieStep(int nComps, int nTypes, int* epsilon, int capacity,
                   int* curState, int* popState, double* clock, long* matrix,
                   long* rowSum, long* totalSum, double rexp, double runi) {
    // find index of a first value in the array, which is larger than threshold
    //      also store the remainder of threshold
    int findLarger(int from, int len, long* vec, double* threshold) {
        double q=*threshold;
        int i=0;
        for(i=from;i<len;i+=1) {
            if(q<vec[i]) {
                *threshold=q;
                return i;
            }
            q-=vec[i];
        }
        return -1;
    }

    // number of rows (and columns) in the transition rate matrix
    int nElems=nComps*nTypes;

    // generate time step
    double dt=rexp/(*totalSum);

    // determine from which district transition will take place
    double q=runi*(*totalSum);
    int fromWhich=findLarger(0,nElems,rowSum,&q);
    int from=fromWhich*nElems;
    int toWhich=findLarger(from,from+nElems,matrix,&q)-from;

    // update states
    curState[fromWhich]-=1;
    popState[fromWhich % nComps]-=1;
    curState[toWhich]+=1;
    popState[toWhich % nComps]+=1;

    // update transition rate matrix
    updateTransitionRateMatrix(fromWhich % nComps,nComps,nTypes,epsilon,
                               capacity,curState,popState,matrix,rowSum,totalSum);
    updateTransitionRateMatrix(toWhich % nComps,nComps,nTypes,epsilon,capacity,
                               curState,popState,matrix,rowSum,totalSum);

    // update time
    *clock+=dt;
}

/*
 * Generates the model's time series (using constant discretization).
 *
 * Input:
 *      dataPoints - number of points in the time series to generate
 *      dt         - discretization step to use
 *      nComps     - number of compartments
 *      nTypes     - number of types
 *      epsilon    - rate of individual transition (depends on type) times 100
 *                   (multiplication is need to avoid using real numbers and
 *                   numerical errors associated with floating point precission
 *                   and to allow for reasonable variability of epsilon values;
 *                   this multiplication also means that time goes 100 times
 *                   faster than without the multiplication)
 *      capacity   - capacity of compartments
 *      curState   - current occupation of the states. Array indexed by state
 *                   id (see info about i).
 *      rng_seed   - random number generator seed
 * Output:
 *      curState   - see input
 *      output     - time series for every district (stored row-wise)
 */
void getSeries(int dataPoints, double dt, int nComps, int nTypes, int* epsilon,
               int capacity, int* curState, int rng_seed, int* output) {
    int i=0;
    int j=0;
    int k=0;
    // initialize GSL random number generator
    gsl_rng_env_setup();
    gsl_rng * rng=gsl_rng_alloc(gsl_rng_taus);
    long seed=(long)rng_seed;
    gsl_rng_set(rng,seed);
    // generate series
    int point=0;
    double clock=0;
    double nClock=dt;
    double rexp=0;
    double runi=0;
    int nElems=nTypes*nComps;
    int* prevState=(int *)malloc(nElems*sizeof(int));
    // Initialize popState
    int* popState=(int *)malloc(nComps*sizeof(int));
    for(i=0;i<nComps;i+=1) {
        popState[i]=0;
        for(j=0;j<nTypes;j+=1) {
            k=i+j*nComps;
            popState[i]+=curState[k];
        }
    }
    // initialize transition rate matrix
    long* matrix=(long *)malloc(nElems*nElems*sizeof(long));
    long* rowSum=(long *)malloc(nElems*sizeof(long));
    long totalSum=0;
    getTransitionRateMatrix(nComps,nTypes,epsilon,capacity,curState,popState,
                            matrix,rowSum,&totalSum);
    // THE MAIN LOOP
    for(point=0;point<dataPoints;point+=1) {// while we need more data points
        while(clock<nClock) {// do until report time hits
            // store previous state (which observed until updated time comes
            memcpy(prevState,curState,nElems*sizeof(int));
            // generate random numbers used by Gillespie method
            rexp=gsl_ran_exponential(rng,1.0);
            runi=gsl_rng_uniform(rng);
            // do a single step using Gillespie method
            gillespieStep(nComps, nTypes, epsilon, capacity, curState, popState,
                          &clock, matrix, rowSum, &totalSum, rexp, runi);
        }
        // copy observed state into output time series
        memcpy(&output[point*nElems],prevState,nElems*sizeof(int));
        nClock+=dt;
    }
    // free allocated memory
    free(prevState);
    free(popState);
    free(matrix);
    free(rowSum);
    // destroy GSL random number generator
    gsl_rng_free(rng);
}

