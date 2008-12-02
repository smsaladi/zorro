
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>


#define EXTERN extern 
#include "utils.h"
#include "trim.h"
#include "tree.h"
#include "weights.h"
#undef EXTERN

#define EXTERN
#include "sample.h"
#undef EXTERN

typedef struct{
  int ki;
  int kj;
  double val;
}keyval;

keyval *sampleNums;
int **sampleFlags;
int do_sample;
int Nsamples;

static int kvcmp(const void *kv1,const void *kv2);

void initSampling(int dummy){
  int g = 0;
  int i,j,r;
  int N = (Nseq*(Nseq-1))/2;
  
  srand(7);

  if(NULL == (sampleNums = (keyval *)(malloc(N*sizeof(keyval))))){
    error("Cannot Allocate Enough Memory\n");
  }
  if(NULL == (sampleFlags = (int **)(malloc(Nseq*sizeof(int *))))){
    error("Cannot Allocate Enough Memory\n");
  }
  for(i=0;i<Nseq;i++){  
    if(NULL == (sampleFlags[i] = (int *)(malloc(Nseq*sizeof(int))))){
      error("Cannot Allocate Enough Memory\n");
    }
    for(j=i+1;j<Nseq;j++){
      sampleFlags[i][j] = 0;
      r = rand();
      sampleNums[g].val =  (((double)(r))/((double)(RAND_MAX)))*((double)(N))*pairWeights[i][j];
      //sampleNums[g].val =  (((double)(r))/((double)(RAND_MAX)));
      sampleNums[g].ki = i;
      sampleNums[g].kj = j;
      g++;
    }
  }

  qsort(sampleNums,g,sizeof(keyval),kvcmp);
  
  if(Nsamples > N){
    Nsamples = N;
  }
  for(i=N-1;i>=N-Nsamples;i--){
    sampleFlags[sampleNums[i].ki][sampleNums[i].kj] = 1;
  }

  for(i=0;i<Nseq;i++){
    fprintf(stderr,"%s: ",leaves[i]->label);
    for(j=0;j<Nseq;j++){
      fprintf(stderr,"%d,%f ",sampleFlags[i][j],pairWeights[i][j]*(double)(N));
    }
    fprintf(stderr,"\n");
  }
  
}


static int kvcmp(const void *v1,const void *v2){
  double r;
  keyval *kv1 = (keyval *)(v1);
  keyval *kv2 = (keyval *)(v2);
  if(kv1->val < kv2->val){
    return -1;
  }
  
  if(kv1->val == kv2->val){
    r = ((double)rand())/((double)RAND_MAX);
    if(r < 0.5){
      return -1;
    }
  }
  
  return 1;
}
