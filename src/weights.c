
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
#include "sample.h"
#undef EXTERN

#define EXTERN
#include "weights.h"
#undef EXTERN

char command[200];
PhyloTree *qtree;
PhyloTree **leaves;
PhyloTree **vlist;
int Nodes;
double **pairWeights;

void reorder_nodes(int dummy);
void initPairWeights(int dummy);
void trace(int base,double l,double WV,PhyloTree *node,PhyloTree *prev);
void trace2(int base,double l,double WV,PhyloTree *node,PhyloTree *prev);
double calcPD( PhyloTree *node );
void calcPDV(PhyloTree *node);

void initWeighting(char *inFile){
  
  char tfile[100];
  int num, errcode;
  double f;
  
  if(Nseq == 1){
    error("Only one sequence in the alignment. No need to run masking!!!!\n");
  }

  if(Nseq == 2){
    int i,j;
    pairWeights = (double **)(malloc(Nseq*sizeof(double *)));
    dists = (double **)(malloc(Nseq*sizeof(double *)));
    for(i=0;i<Nseq;i++){
      pairWeights[i] = (double *)(malloc(Nseq*sizeof(double)));
      dists[i] = (double *)(malloc(Nseq*sizeof(double)));
      for(j=0;j<Nseq;j++){
	if(j == i){
	  pairWeights[i][j] = 0.0;
	  dists[i][j] = 0.0;
	}
	else{
	  pairWeights[i][j] = 1.0;
	  dists[i][j] = 1.0;
	}
      }
    }
    return;
  }
  if(uguide == 0){
    sprintf(tfile,"%s.zorro.tre",inFile);
    sprintf(command,"%s %s > %s 2>/dev/null",treeprog,inFile,tfile);
    errcode=system(command);
    qtree = readPhyloTree(tfile);
    sprintf(command,"rm -rf %s",tfile);
    errcode=system(command);
  }
  else{
    qtree = readPhyloTree(guidetree);
  }

  
  
  if(qtree == NULL) {
    error("Cannot create tree from alignment \n Check %s installation\n",treeprog);
  }
  
  makeBinary(qtree);
  if(qtree->numChild == 2){
    f = qtree->child[0]->ewt + qtree->child[1]->ewt;
    qtree->child[0]->ewt = f;
    qtree->child[1]->ewt = f;
  }
  
  leaves = getAllLeaves(qtree,&num);
  
  if(num != Nseq){
    error("Tree inconsistent with alignment\n Check %s installation\n",treeprog);
  }
  // vlist in postorder
  vlist = getAllNodes(qtree,&Nodes);
 
  reorder_nodes(0);
  initPairWeights(0);
}


#ifdef RATIONALE2

void initPairWeights(int dummy){
  int i,j,k;
  PhyloTree *no;
  PhyloTree *bro,*par;
  double prod,sum,x;
  for(k=0;k<Nodes;k++){
    if(vlist[k]->ewt<=1e-8){
      vlist[k]->ewt = 1e-8;
    }
    vlist[k]->ewt += 1e-8;
  }
  // vlist in postorder
  for(k=0;k<Nodes-1;k++){
    no = vlist[k];
    if(no->numChild == 0){
      no->w = log(1.0);
      no->W = log(no->ewt);
    }
    else{
      no->w = -FLT_MAX;
      x = 0.0;
      if(k == Nodes-2){
	bro = no->parent->child[0];
      }

      for(i=0;i<no->numChild;i++){	
	prod = no->child[i]->w;
	x += no->child[i]->W;
	for(j=0;j<no->numChild;j++){
	  if(i!=j){
	    prod = prod + no->child[j]->W;   
	  }
	}
	if(k == Nodes-2){
	  prod = prod + bro->w;
	}
	addLogProb(prod,no->w);
      }    
     
      if(k == Nodes-2){
	prod = bro->w;
	x += bro->W;
	for(j=0;j<no->numChild;j++){
	  prod = prod + no->child[j]->W;  
	}
	addLogProb(prod,no->w);
      }
      
      no->W = log(no->ewt) + no->w;
      addLogProb(x,no->W);
    }
  }

  no=vlist[Nodes-2];
  no->v = -FLT_MAX;
  no->V = 0;
  
  for(k=Nodes-3;k>=0;k--){
    no = vlist[k];
    par = no->parent;
    prod = 0.0;
    sum = -FLT_MAX;
    // special case for root's other child
    if(no != vlist[Nodes-1]->child[0]){     
      for(i=0;i<par->numChild;i++){
	if(no == par->child[i]){
	  continue;
	}
	bro = par->child[i];
	prod += bro->W;
	addLogProb(bro->w-bro->W,sum);
      }
    
      if(par == vlist[Nodes-2]){
	bro = vlist[Nodes-1]->child[0];
	prod += bro->W;
	addLogProb(bro->w-bro->W,sum);
      }
      else{
	prod += par->V;
	addLogProb(par->v-par->V,sum);
      }
    }
    else{
      par = vlist[Nodes-2];
      for(i=0;i<par->numChild;i++){
	bro = par->child[i];
	prod += bro->W;
	addLogProb(bro->w-bro->W,sum);
      }
    }
    no->v = sum + prod;
    no->V += log(no->ewt) + no->v;
    addLogProb(prod,no->V);
  }
  
  pairWeights = (double **)(malloc(Nseq*sizeof(double *)));
  dists = (double **)(malloc(Nseq*sizeof(double *)));
  for(i=0;i<Nseq;i++){
    pairWeights[i] = (double *)(malloc(Nseq*sizeof(double)));
    dists[i] = (double *)(malloc(Nseq*sizeof(double)));
    for(j=0;j<Nseq;j++){
      pairWeights[i][j] = 0.0;
      dists[i][j] = 0.0;
    }
  }

  

  for(i=0;i<Nseq;i++){
    if(leaves[i] != vlist[Nodes-1]->child[0]){
      trace(i,leaves[i]->ewt,0.0,leaves[i]->parent,leaves[i]);
    }
    else{
      trace(i,leaves[i]->ewt,0.0,vlist[Nodes-2],leaves[i]);
    }
  }
  
  sum = 0.0;
  for(i=0;i<Nseq;i++){
    for(j=0;j<Nseq;j++){
      sum += pairWeights[i][j];
    }
  }
  sum = sum / 2;

  for(i=0;i<Nseq;i++){
    if(!do_sample){
      if(verbose)
	fprintf(stderr,"%s: ",leaves[i]->label);
    }
    for(j=0;j<Nseq;j++){
      pairWeights[i][j] /= sum;
      if(!do_sample){
	if(verbose)
	  fprintf(stderr,"%f ",pairWeights[i][j]);
      }
    }
    if(!do_sample){
      if(verbose)
	fprintf(stderr,"\n");
    }

}


void trace(int base,double l,double WV,PhyloTree *node,PhyloTree *prev){
  int i;
  double pi;
  

  if(node->numChild == 0){
    pairWeights[base][node->index] = l * l * l *exp(WV);
    dists[base][node->index] = l;
    return;
  }
  
  
  /* Special cases == children of root */
  if(node == vlist[Nodes-2]){
    
    pi = 0.0;
    for(i=0;i<node->numChild;i++){
      if(node->child[i] != prev){
	pi += node->child[i]->W;
      }
    }
    
    if(prev != vlist[Nodes-1]->child[0]) {
      pi += vlist[Nodes-1]->child[0]->W; 
    }
  
    for(i=0;i<node->numChild;i++){
      if(node->child[i] != prev){
	trace(base,l+node->child[i]->ewt,WV+pi-node->child[i]->W,node->child[i],node);
      }
    }
    if(prev != vlist[Nodes-1]->child[0]) {
      trace(base,l+node->ewt,WV+pi-vlist[Nodes-1]->child[0]->W,vlist[Nodes-1]->child[0],node);
    }  

    return;
  }
  
  if(node == vlist[Nodes-1]->child[0]){
    pi = 0.0;
    for(i=0;i<node->numChild;i++){
      if(node->child[i] != prev){
	pi += node->child[i]->W;
      }
    }
    
    if(prev != vlist[Nodes-2]) {
      pi += node->V; 
    }
  

    for(i=0;i<node->numChild;i++){
      if(node->child[i] != prev){
	trace(base,l+node->child[i]->ewt,WV+pi-node->child[i]->W,node->child[i],node);
      }
    }

    if(prev != vlist[Nodes-2]) {
      trace(base,l+node->ewt,WV+pi-node->V,vlist[Nodes-2],node);
    }
    

    return;
  }

  /*All other cases*/
  pi = 0.0;
  for(i=0;i<node->numChild;i++){
    if(node->child[i] != prev){
      pi += node->child[i]->W;
    }
  }
  
  if(prev != node->parent) {
    pi += node->V; 
  }
  
  for(i=0;i<node->numChild;i++){
      if(node->child[i] != prev){
	trace(base,l+node->child[i]->ewt,WV+pi-node->child[i]->W,node->child[i],node);
      }
    }

  if(prev != node->parent) {
    trace(base,l+node->ewt,WV+pi-node->V,node->parent,node);
  }
  
  
}

#else 
#ifdef RATIONALE1
void initPairWeights(int dummy){
  int i,j;
  double sum;

  pairWeights = (double **)(malloc(Nseq*sizeof(double *)));
  dists = (double **)(malloc(Nseq*sizeof(double *)));
  for(i=0;i<Nseq;i++){
    pairWeights[i] = (double *)(malloc(Nseq*sizeof(double)));
    dists[i] = (double *)(malloc(Nseq*sizeof(double)));
    for(j=0;j<Nseq;j++){
      pairWeights[i][j] = 0.0;
      dists[i][j] = 0.0;
    }
  }

  for(i=0;i<Nseq;i++){
    if(leaves[i] != vlist[Nodes-1]->child[0]){
      trace(i,leaves[i]->ewt,0.0,leaves[i]->parent,leaves[i]);
    }
    else{
      trace(i,leaves[i]->ewt,0.0,vlist[Nodes-2],leaves[i]);
    }
  }
  
  sum = 0.0;
  for(i=0;i<Nseq;i++){
    for(j=0;j<Nseq;j++){
      sum += pairWeights[i][j];
    }
  }
  sum = sum / 2;

  for(i=0;i<Nseq;i++){
    if(verbose)
      fprintf(stderr,"%s: ",leaves[i]->label);
    for(j=0;j<Nseq;j++){
      pairWeights[i][j] /= sum;
      if(verbose)
	fprintf(stderr,"%f ",pairWeights[i][j]);
    }
    if(verbose)
      fprintf(stderr,"\n");
  }


}

void trace(int base,double l,double WV,PhyloTree *node,PhyloTree *prev){
  int i;
  

  if(node->numChild == 0){
    pairWeights[base][node->index] = l * l *exp(WV);
    dists[base][node->index] = l;
    return;
  }
  
  
  /* Special cases == children of root */
  if(node == vlist[Nodes-2]){
     
  
    for(i=0;i<node->numChild;i++){
      if(node->child[i] != prev){
	trace(base,l+node->child[i]->ewt,WV-log(node->numChild),node->child[i],node);
      }
    }
    if(prev != vlist[Nodes-1]->child[0]) {
      trace(base,l+node->ewt,WV-log(node->numChild),vlist[Nodes-1]->child[0],node);
    }  

    return;
  }
  
  if(node == vlist[Nodes-1]->child[0]){
    
    for(i=0;i<node->numChild;i++){
      if(node->child[i] != prev){
	trace(base,l+node->child[i]->ewt,WV-log(node->numChild),node->child[i],node);
      }
    }

    if(prev != vlist[Nodes-2]) {
      trace(base,l+node->ewt,WV-log(node->numChild),vlist[Nodes-2],node);
    }
    

    return;
  }

  /*All other cases*/
  
  
  for(i=0;i<node->numChild;i++){
    if(node->child[i] != prev){
      trace(base,l+node->child[i]->ewt,WV-log(node->numChild),node->child[i],node);
    }
  }
  
  if(prev != node->parent) {
    trace(base,l+node->ewt,WV-log(node->numChild),node->parent,node);
  }
  
  
}

#else 
#ifdef ALGO1

void initPairWeights(int dummy){
  int i,j;
  double sum;
  
  calcDescLeaves(qtree);
  
  pairWeights = (double **)(malloc(Nseq*sizeof(double *)));
  dists = (double **)(malloc(Nseq*sizeof(double *)));
  for(i=0;i<Nseq;i++){
    pairWeights[i] = (double *)(malloc(Nseq*sizeof(double)));
    dists[i] = (double *)(malloc(Nseq*sizeof(double)));
    for(j=0;j<Nseq;j++){
      pairWeights[i][j] = 0.0;
      dists[i][j] = 0.0;
    }
  }

  for(i=0;i<Nseq;i++){
    double d = (double)(Nseq-1);
    if(leaves[i] != vlist[Nodes-1]->child[0]){
      trace(i,leaves[i]->ewt,leaves[i]->ewt/d,leaves[i]->parent,leaves[i]);
    }
    else{
      trace(i,leaves[i]->ewt,leaves[i]->ewt/d,vlist[Nodes-2],leaves[i]);
    }
  }
  
  sum = 0.0;
  for(i=0;i<Nseq;i++){
    for(j=0;j<Nseq;j++){
      sum += pairWeights[i][j];
    }
  }
  sum = sum / 2;

  for(i=0;i<Nseq;i++){
    if(!do_sample){
      if(verbose)
	fprintf(stderr,"%s: ",leaves[i]->label);
    }
    for(j=0;j<Nseq;j++){
      pairWeights[i][j] /= sum;
      if(!do_sample){
	if(verbose)
	  fprintf(stderr,"%f ",pairWeights[i][j]);
      }
      
      if(!do_sample){
	if(verbose)
	  fprintf(stderr,"\n");
      }
    }
    
    
  }
}

void trace(int base,double l,double WV,PhyloTree *node,PhyloTree *prev){
  int i;
  double Np;
  PhyloTree *up;

  if(node->numChild == 0){
    pairWeights[base][node->index] = WV;
    dists[base][node->index] = l;
    return;
  }
 
  /* Special cases == children of root */
  if(node == vlist[Nodes-2]){
    up = vlist[Nodes-1]->child[0];
  }
  else if(node == vlist[Nodes-1]->child[0]) {
    up = vlist[Nodes-2];
  }else{/*All other cases*/
    up = node->parent;
  }
  /* Recursion*/
  for(i=0;i<node->numChild;i++){
    Np = ((double)(node->child[i]->DescLeaves))*((double)(Nseq-node->child[i]->DescLeaves));
    if(node->child[i] != prev){
      trace(base,l+node->child[i]->ewt,WV+(node->child[i]->ewt/Np),node->child[i],node);
    }
  }
  
  if(prev != up) {
    Np = ((double)(node->DescLeaves))*((double)(Nseq-node->DescLeaves));
    trace(base,l+node->ewt,WV+(node->ewt/Np),up,node);
  }
  
  
}

#else

void initPairWeights(int dummy){
  int i,j,k;
  double sum;
  
  for(k=0;k<Nodes;k++){
    if(vlist[k]->ewt<=1e-8){
      vlist[k]->ewt = 1e-8;
    }
    vlist[k]->ewt += 1e-8;
  }

  calcDescLeaves(qtree);
  calcPD(qtree);
  calcPDV(vlist[Nodes-2]);
  calcPDV(vlist[Nodes-1]->child[0]);
  
  pairWeights = (double **)(malloc(Nseq*sizeof(double *)));
  dists = (double **)(malloc(Nseq*sizeof(double *)));
  for(i=0;i<Nseq;i++){
    pairWeights[i] = (double *)(malloc(Nseq*sizeof(double)));
    dists[i] = (double *)(malloc(Nseq*sizeof(double)));
    for(j=0;j<Nseq;j++){
      pairWeights[i][j] = 0.0;
      dists[i][j] = 0.0;
    }
  }

  for(i=0;i<Nseq;i++){
    if(leaves[i] != vlist[Nodes-1]->child[0]){
      trace2(i,log(leaves[i]->ewt),log(1.0),leaves[i]->parent,leaves[i]);
    }
    else{
      trace2(i,log(leaves[i]->ewt),log(1.0),vlist[Nodes-2],leaves[i]);
    }
  }
  
  sum = 0.0;
  for(i=0;i<Nseq;i++){
    for(j=0;j<Nseq;j++){
      sum += pairWeights[i][j];
    }
  }
  sum = sum / 2;

  for(i=0;i<Nseq;i++){
    if(!do_sample){
      if(verbose)
	fprintf(stderr,"%s: ",leaves[i]->label);
    }
    for(j=0;j<Nseq;j++){
      if(!do_sample){
	if(verbose)
	  fprintf(stderr,"%f ",pairWeights[i][j]);
      }
    
      pairWeights[i][j] /= sum;
    }
    if(!do_sample){
      if(verbose)
	fprintf(stderr,"\n");
    }
  }


}

void trace2(int base,double l,double WV,PhyloTree *node,PhyloTree *prev){
  int i;
  double r1,r2;
  PhyloTree *up;

  if(node->numChild == 0){
    pairWeights[base][node->index] = exp(l/2);
    return;
  }
 
  /* Special cases == children of root */
  if(node == vlist[Nodes-2]){
    up = vlist[Nodes-1]->child[0];
  }
  else if(node == vlist[Nodes-1]->child[0]) {
    up = vlist[Nodes-2];
  }else{/*All other cases*/
    up = node->parent;
  }



  /* Recursion*/  

  if(prev != up) {
    // Go up
    r1 = -FLT_MAX;
    for(i=0;i<node->numChild;i++){
      addLogProb(log(node->child[i]->w),r1);
    }
    r1 = log(prev->w)-r1;
    r1 = WV+r1;
    
    r2 = log(node->v);
    for(i=0;i<node->numChild;i++){
      if(prev != node->child[i])
	addLogProb(log(node->child[i]->w),r2);
    }
  
    r2 = log(node->v)+ l - r2; 
    addLogProb(r1+log(node->ewt),r2);
    trace2(base,r2,r1,up,node);
    // Go down
    r1 = log(prev->w);
    addLogProb(log(node->v),r1);
    r1 = log(prev->w)-r1;
    r1 = WV+r1;
    
    for(i=0;i<node->numChild;i++){
      if(prev != node->child[i]){
	//r2 = (node->child[i]->w)/(node->child[i]->w+node->v);
	//r2 = r2*l + r1*node->child[i]->ewt;
	r2 = log(node->v);
	addLogProb(log(node->child[i]->w),r2);
	r2 = log(node->child[i]->w) - r2;
	r2 = r2+l;
	addLogProb(r1+log(node->child[i]->ewt),r2);
	
	trace2(base,r2,r1,node->child[i],node);
      }
    }
  }
  else{
    // Can only go down
    for(i=0;i<node->numChild;i++){
      int j;
      r1 = log(node->v);
      for(j=0;j<node->numChild;j++){
	if(i != j){
	  addLogProb(log(node->child[j]->w),r1);
	}
      }
      
      r1 = (log(node->v)-r1)+WV;
      
      r2 = -FLT_MAX;
      for(j=0;j<node->numChild;j++){
	addLogProb(log(node->child[j]->w),r2);
      }
      r2 = (log(node->child[i]->w)-r2);
      r2 = r2+l;
      addLogProb(r1+log(node->child[i]->ewt),r2);
      
      trace2(base,r2,r1,node->child[i],node);   
    }
  }
}
 
  
 



double calcPD( PhyloTree *node ){
  int i;
  if( node->numChild == 0 ){
    node->W = 0.0;
    node->w = node->ewt;
    return node->w;
  }
  node->W = 0.0;
  for( i = 0; i < node->numChild; i++ ){
    node->W += calcPD(node->child[i]);    
  }
  node->w = node->W + node->ewt;
  return node->w;
}

void calcPDV(PhyloTree *node){
  PhyloTree *up;
  int i;
  /* Special cases == children of root */
  if(node == vlist[Nodes-2]){
    up = vlist[Nodes-1]->child[0];
    node->V = up->W;
    node->v = up->w;
    for( i = 0; i < node->numChild; i++ )
      calcPDV(node->child[i]);
    return;
  }
  else if(node == vlist[Nodes-1]->child[0]) {
    up = vlist[Nodes-2];
    node->V = up->W;
    node->v = up->w;
    for( i = 0; i < node->numChild; i++ )
      calcPDV(node->child[i]);
    return;
  }

  /*All other cases*/
  up = node->parent;
  node->V = up->v;
  for( i = 0; i < up->numChild; i++ )
    if(up->child[i] != node)
      node->V += up->child[i]->w;
  
  node->v = node->ewt + node->V;
  for( i = 0; i < node->numChild; i++ )
      calcPDV(node->child[i]);
  return;
}

#endif
#endif
#endif




void reorder_nodes(int dummy){
  int i,j;
  PhyloTree *temp;
  for(i=0;i<Nseq;i++){
    for(j=i;j<Nseq;j++){
      if(strcmp(names[i],leaves[j]->label)==0){
	leaves[j]->index = i;
	break;
      }
    }
    if(j==Nseq){
      error("Name %s not found in the Tree Leaves\n",names[i]);
    }
    
    if(j!=i){
      temp = leaves[j];
      leaves[j] = leaves[i];
      leaves[i] = temp;
    }
    //printf("%d %d %s\n",i,j,leaves[i]->label);
  }
}
