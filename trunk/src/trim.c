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

#ifdef WIN32
int strcasecmp(const char *s1, const char *s2){ return _strcmpi(s1, s2); }
#endif

int main(int argc,char **argv){ 
  int i,t;
  if(argc < 2){
    printf("File not specified\n");
    printHelp(0);
    exit(EXIT_FAILURE);
  }

  if(0 == strcasecmp(argv[1],"-help") || 0 == strcasecmp(argv[1],"-h")){
      printHelp(0);
      exit(EXIT_FAILURE);
  }

  // User provided options

  do_sample = 0;
  uguide = 0;
  Nsamples = 10*Nseq;
  strcpy(treeprog,"FastTree");
  ignoreGaps = 0;
  WEIGHTING = 1;
  verbose = 0;
  // Hard coded options
  JTT = 0;
  PMB = 0;
  PAM = 1;
  MATRICES = 0;
  ADD = 1;

  for(i=1;i<argc-1;i++){
    if(0 == strcasecmp(argv[i],"-help") || 0 == strcasecmp(argv[i],"-h")){
      printHelp(0);
    }
    
    if(0 == strcasecmp(argv[i],"-verbose") || 0 == strcasecmp(argv[i],"-v")){
      verbose = 1;
    }
    
    if(0 == strcasecmp(argv[1],"-sample")){
      do_sample = 1;
      if(verbose)
	fprintf(stderr,"Sampling\n");
    }

    if(0 == strcasecmp(argv[i],"-nosample")){
      do_sample = 0;
      if(verbose)
	fprintf(stderr,"All Pairs Used No Sampling... might take longer time\n");
    }
    
    if(0 == strcasecmp(argv[1],"-ignoregaps")){
      ignoreGaps = 1;
      if(verbose)
	fprintf(stderr,"Ignoring gap-pairs\n");
    }

    if(0 == strcasecmp(argv[1],"-noweighting")){
      WEIGHTING = 0;
      if(verbose)
	fprintf(stderr,"Using sum of pairs scheme instead of weighted sum of pairs\n");
    }



    if(0 == strcasecmp(argv[i],"-Nsample")){
      do_sample = 1;
      if(verbose)
	fprintf(stderr,"User provided Nsamples\n");
      if(i == argc-2){
	error("Number of pairs to be sampled not provided\n");
      }
      else{
	i++;
	sscanf(argv[i],"%d%n",&Nsamples,&t);
	if(t==0){
	  error("Nsamples provided not an integer\n");
	}
	else{
	  if(verbose)
	    fprintf(stderr,"Sampling %d pairs\n",Nsamples);
	}
      }
    }

    if(0 == strcasecmp(argv[i],"-treeprog")){
      if(verbose)
	fprintf(stderr,"User provided guide Tree Inference Program\n");
      if(i == argc-2){
	error("Tree inference program not provided\n");
      }
      else{
	i++;
	strcpy(treeprog,argv[i]);
      }
    }

    if(0 == strcasecmp(argv[i],"-guide")){
      uguide = 1;
      if(verbose)
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
    if(Nsamples < 10*Nseq){
      Nsamples = 10*Nseq; 
    }
    initSampling(Nsamples);
  }

  initHMM(alen);
  calc_posterior(alen);

  return 1;
  
}

void printHelp(int dummy){
  fprintf(stdout,"Usage: zorro [options] filename > outputfile\n\n");
      fprintf(stdout,"ZORRO Options \n\n");
      fprintf(stdout,"-sample          : Sampling pairs to calculate alignment reliabilty [Off By Default]\n");
      fprintf(stdout,"-nosample        : No Sampling i.e. using every pair to calculate alignment reliabilty [On By Default]\n");
      fprintf(stdout,"-noweighting     : Using sum of pairs instead of weighted sum of pairs to calculate column confidence [Off By Default]\n");
      fprintf(stdout,"-ignoregaps      : Ignore pair-gaps in columns when calculating column confidences [Off By Default]\n");
      fprintf(stdout,"-Nsample NUMBER  : Tells ZORRO to sample #NUMBER pairs when sampling, automatically turns on -sample option [Samples 10*Nseq sequences By Default]\n");
      fprintf(stdout,"-treeprog PROGRAM: Tells ZORRO to use PROGRAM instead of the default FastTree to create guide tree [FastTree By Default]\n");
      fprintf(stdout,"-guide treefile  : User provided guide tree\n");     
      fprintf(stdout,"\n");
      exit(0);
}
