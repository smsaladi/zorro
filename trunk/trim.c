#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>


#define EXTERN
#include "trim.h"
#undef EXTERN

#define EXTERN extern
#include "utils.h"
#include "hmm.h"
#include "tree.h"
#include "weights.h"
#include "sample.h"
#include "matrices.h"
#undef EXTERN

void printHelp(int dummy);

int main(int argc,char **argv){ 
  int i;
  if(argc < 2){
    printf("File not specified\n");
    printHelp(0);
    exit(EXIT_FAILURE);
  }


  if(0 == strcasecmp(argv[1],"-help") || 0 == strcasecmp(argv[1],"-h")){
      printHelp(0);
  }

  // User provided options

  do_sample = 0;
  uguide = 0;

  // Hard coded options
  JTT = 0;
  PMB = 0;
  PAM = 1;
  MATRICES = 0;

  for(i=1;i<argc-1;i++){
    if(0 == strcasecmp(argv[i],"-help") || 0 == strcasecmp(argv[i],"-h")){
      printHelp(0);
    }
    
    if(0 == strcmp(argv[1],"-sample")){
      do_sample = 1;
      fprintf(stderr,"Sampling\n");
    }

    if(0 == strcasecmp(argv[i],"-nosample")){
      do_sample = 0;
      fprintf(stderr,"All Pairs Used No Sampling... might take longer time\n");
    }

    if(0 == strcasecmp(argv[i],"-guide")){
      uguide = 1;
      fprintf(stderr,"User provided guide tree\n");
      if(i == argc-2){
	error("Guide tree file-name not provided\n");
      }
      else{
	i++;
	strcpy(guidetree,argv[i]);
      }
    }
  }

  readSeq(argv[argc-1]);
  initWeighting(argv[argc-1]);
  if(do_sample){
    Nsamples = 4*Nseq; 
    initSampling(Nsamples);
  }

  initHMM(alen);
  calc_posterior(alen);

  return 1;
  
}

void printHelp(int dummy){
  fprintf(stdout,"Usage: probmask [options] filename\n\n");
      fprintf(stdout,"ZORRO Options \n\n");
      fprintf(stdout,"-sample         : Sampling pairs to calculate alignment reliabilty\n");
      fprintf(stdout,"-nosample       : No Sampling i.e. using every pair to calculate alignment reliabilty\n");
      fprintf(stdout,"-guide treefile : User provided guide tree\n");     
      fprintf(stdout,"\n");
      exit(0);
}
