/*
Routines partially adapted from AVID source code:
N. Bray, I. Dubchak, and L. Pachter, AVID: A Global Alignment Program, Genome Research 13, 97-102 (2003).
*/


#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>

#define EXTERN extern
#include "trim.h"
#include "utils.h"
#undef EXTERN

#include "tree.h"


#define MAX_LEN 30

void deletePhyloTree( PhyloTree *tree )
{
  int i;
  for(i = 0; i < tree->numChild; i++ ){
    deletePhyloTree( tree->child[i] );
  }
  if( tree->numChild > 0 ){
    free(tree->child);
  }
  free(tree);
  
}


void deleteNode( PhyloTree *tree )
{
  
  if( tree->numChild > 0 ){
    free(tree->child);
  }
  free(tree);

}


PhyloTree *copyPhyloTree( PhyloTree *tree )
{

  PhyloTree *toRet = (PhyloTree *)(malloc(sizeof(PhyloTree)));
  int i;;
  if( tree->label ){
    toRet->label = (char *)(malloc(sizeof(char)*(strlen(tree->label) + 1)));
    strcpy( toRet->label, tree->label );
  }
  else {
    toRet->label = NULL;
  }

  toRet->parent = tree->parent;
  toRet->numChild = tree->numChild;
  toRet->prof = tree->prof;
  toRet->ewt = tree->ewt;

  if( toRet->numChild > 0 ){
    toRet->child = (PhyloTree **)(malloc(sizeof(PhyloTree *)*toRet->numChild));
    for(i = 0; i < toRet->numChild; i++ ){
      toRet->child[i] = copyPhyloTree( tree->child[i] );
      toRet->child[i]->parent = toRet;
    }
  }

  return toRet;
  
}

PhyloTree *readPhyloTree( char *fileName )
{

  int slen;
  char *string;
  FILE *in;
  PhyloTree *toRet;

  in = fopen( fileName, "r" );
  if( !in ){
    error("Cannot open Tree File %s for reading\n",fileName);
  }
  
  for( slen = 0; !feof( in ); ){
    char c = getc( in );

    if( !isspace( c ) && c != EOF )
      slen++;
  }
  string = (char *)(malloc((slen + 1)*sizeof(char)));



  rewind( in );
  for( slen = 0; !feof( in ); ){
    char c = getc( in );

    if( !isspace( c ) && c != EOF )
      string[slen++] = c;
  }
  string[slen] = '\0';

  if( string[0] != '(' ){
    free(string);
    fclose( in );
    return NULL;
  }
  
  
  toRet = stringToPhyloTree( string + 1);

  free(string);
  fclose( in );

  if(toRet == NULL){
    fprintf(stderr,"Warning : Empty Tree from Non Empty File %s\n",fileName);
  }

  return toRet;

}

PhyloTree *stringToPhyloTree( char *string )
{
  int i,j,depth,N;
  PhyloTree *toRet = (PhyloTree *)(malloc(sizeof(PhyloTree)));
  toRet->numChild = 1;
  toRet->label = NULL;
  toRet->index = -1;
  toRet->parent = NULL;
  toRet->prof = NULL;
  toRet->ewt = 0.0;

  for( depth = 0, i = 0; depth >= 0; i++ ){
    if( string[i] == '\0' )
      return NULL;
    if( string[i] == '(' )
      depth++;
    if( string[i] == ')' )
      depth--;
    if( string[i] == ',' && depth == 0 )
      toRet->numChild++;
  }

  //printf("Children = %d \n",toRet->numChild);

  
  toRet->child = (PhyloTree **)(malloc(toRet->numChild*sizeof(PhyloTree *)));
  for( j = 0; j < toRet->numChild; j++ ){
    toRet->child[j] = NULL;
  }
  toRet->parent = NULL;

  N = toRet->numChild;
  toRet->numChild = 0;
  for( i = 0; string[i] != ')'; ){
    if( string[i] == '(' ){
      toRet->child[toRet->numChild] = stringToPhyloTree( string + i + 1 );
      if( !toRet->child[toRet->numChild] ){
	fprintf(stderr,"Warning : Probable Inconsistency in reading number of Sub-Trees %d %d\n",N,toRet->numChild);
	for( j = 0; j < toRet->numChild; j++ ){
	  free(toRet->child[j]);
	}
	free(toRet);	
	return NULL;
      }
      toRet->child[toRet->numChild]->parent = toRet;
      if( toRet->child[toRet->numChild] == NULL )
	return NULL;
      i++;
      for( depth = 1; depth > 0; i++ ){
	if( string[i] == '(' )
	  depth++;
	if( string[i] == ')' )
	  depth--;
      }
      // Read bootstrap values
      if( string[i] != ':' ){
	int temp = i; 
	while( string[i] == '.' || isdigit(string[i]) ){
	  i++;
	}
	if(string[i] != ':'){
	  fprintf(stderr,"Warning : Probable Inconsistency, expected \":\" not found\n");
	  for( j = 0; j < toRet->numChild; j++ ){
	    free(toRet->child[j]);
	  }
	  free(toRet);
	  return NULL;
	}
	toRet->child[toRet->numChild]->bootstrap =  (double)atof(string + temp);
      }
      else{
	toRet->child[toRet->numChild]->bootstrap = -1;
      }
      sscanf(string + i + 1, "%lf" , &toRet->child[toRet->numChild]->ewt);
      toRet->numChild++;
      while( string[i] != ',' && string[i] != ')' )
	i++;
    }
    else {
      int temp = i;

      while( string[i] != ':' && string[i] != ',' && string[i] != ')' )
	i++;
      if( string[i] != ':' || temp == i )
	return NULL;
      toRet->child[toRet->numChild] = (PhyloTree *)(malloc(sizeof(PhyloTree)));
      toRet->child[toRet->numChild]->numChild = 0;
      toRet->child[toRet->numChild]->ewt = 0.0;
      toRet->child[toRet->numChild]->parent = toRet;
      toRet->child[toRet->numChild]->child = NULL;
      toRet->child[toRet->numChild]->prof = NULL;
      toRet->child[toRet->numChild]->index = -1;
      toRet->child[toRet->numChild]->label = (char *)(malloc(sizeof(char)*(i - temp + 1)));
      i = temp;
      while( string[i] != ':' ){
	toRet->child[toRet->numChild]->label[i - temp] = string[i];
	i++;
      }
	
      toRet->child[toRet->numChild]->label[i - temp] = '\0';
      sscanf(string + i + 1, "%lf" , &toRet->child[toRet->numChild]->ewt);
      toRet->numChild++;
      while( string[i] != ',' && string[i] != ')' )
	i++;
    }
    if( string[i] == ',' )
      i++;
  }
  
  //printf("Children = %d \n",toRet->numChild);

  return toRet;

}




int isBinaryPhyloTree( PhyloTree *tree )
{

  if( !tree )
    return FALSE;

  if( tree->numChild == 0 )
    return TRUE;

  if( tree->numChild != 2 )
    return FALSE;

  return isBinaryPhyloTree( tree->child[0] ) &&
    isBinaryPhyloTree( tree->child[1] );

}

void makeBinary(PhyloTree *node){
  PhyloTree **tlist;
  int i;
  
 
  
  if( node->numChild <= 2 ){
    if(node->numChild == 0){
      return;
    }
    if(node->numChild == 1){
      fprintf(stderr,"Warning : Not yet implemented module for making node with one child binary\n");   
    }
    for(i=0;i<node->numChild;i++){
      makeBinary(node->child[i]);
    }
    /* Make sure that the last child of the root is not a leaf*/
    if(node->parent == NULL){
      if(node->child[1]->numChild == 0){
	PhyloTree *tmp;
	tmp = node->child[0];
	node->child[0] = node->child[1];
	node->child[1] = tmp;
      }
    }
    return;
  }
  tlist = node->child;
  node->numChild -= 1;
  i = node->numChild-1;
  node->child = (PhyloTree **)(malloc(node->numChild*sizeof(PhyloTree *)));
  node->child[i] = (PhyloTree *)(malloc(sizeof(PhyloTree)));  
  node->child[i]->ewt = 1e-8;
  node->child[i]->label = NULL;
  node->child[i]->index = -1; 
  node->child[i]->parent = node;
  node->child[i]->numChild = 2;
  node->child[i]->child = (PhyloTree **)(malloc(2*sizeof(PhyloTree *)));
  node->child[i]->child[0] = tlist[i];
  node->child[i]->child[1] = tlist[i+1];
  node->child[i]->child[0]->parent =  node->child[i];
  node->child[i]->child[1]->parent =  node->child[i];
  
  for(i=0;i<node->numChild-1;i++){
    node->child[i] = tlist[i];
  }  
  free(tlist);
  for(i=0;i<node->numChild;i++){
    makeBinary(node->child[i]);
  }

}


int countPhyloTreeNodes( PhyloTree *tree )
{
  int i;
  int num = 0;
  if( !tree )
    return 0;

  if( tree->numChild == 0 )
    return 1;
  
  for( i = 0; i < tree->numChild; i++ )
    num += countPhyloTreeNodes( tree->child[i] );

  return num;

}


PhyloTree **getAllNodes( PhyloTree *tree , int *num){
  PhyloTree **nodes;
  int N;

  N = 2*countPhyloTreeNodes( tree );
  *num = N;
  
  nodes = (PhyloTree **)(malloc(N*sizeof(PhyloTree *)));
  *num = 0;
  get_nodes( tree, nodes, num );
  
  return nodes;

}

void get_nodes( PhyloTree *tree, PhyloTree **nodes, int *num )
{
  int i;
  if( !tree )
    return;

 
  for( i = 0; i < tree->numChild; i++ )
    get_nodes( tree->child[i], nodes, num );

  nodes[*num] = tree;
  *num = *num + 1;
}






PhyloTree **getAllLeaves( PhyloTree *tree , int *num)
{

  PhyloTree **leaves;
  int N;

  N = countPhyloTreeNodes( tree );
  *num = N;
  
  leaves = (PhyloTree **)(malloc(N*sizeof(PhyloTree *)));
  *num = 0;
  get_leaves( tree, leaves, num );

  return leaves;

}


void get_leaves( PhyloTree *tree, PhyloTree **leaves, int *num )
{
  int i;
  if( !tree )
    return;

  if( tree->numChild == 0 ){
    leaves[*num] = tree;
    *num = *num + 1;
  }

  for( i = 0; i < tree->numChild; i++ )
    get_leaves( tree->child[i], leaves, num );

}

void calcDescLeaves(PhyloTree *tree){
  int i;
  if(tree->numChild == 0){
    tree->DescLeaves = 1;
    return;
  }
  tree->DescLeaves = 0;
  for(i=0;i<tree->numChild;i++){
    calcDescLeaves(tree->child[i]);
    tree->DescLeaves += tree->child[i]->DescLeaves;
  }

}


int getPhyloTreeIndices( PhyloTree *tree, char **labels, int N)
{
  int i;
  tree->prof = (profile *)(malloc(sizeof(profile)));
  tree->prof->score = NULL;
  tree->prof->index = -1;

  if( tree->numChild == 0 ){
    for( i = 0; i < N; i++ ){
      if( strncmp( tree->label, labels[i], MAX_LABEL_LEN ) == 0){
	tree->prof->index = -1;
      }
    }
    if( tree->prof->index == -1 )
      return FALSE;
  }
    
  for( i = 0; i < tree->numChild; i++ )
    if( !getPhyloTreeIndices( tree->child[i], labels , N) )
      return FALSE;

  return TRUE;

}


