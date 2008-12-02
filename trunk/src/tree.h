#ifndef PHYLOTREE
#define PHYLOTREE

#define FALSE 0
#define TRUE 1


struct PhyloTree;


typedef struct{
  struct PhyloTree *node;
  int index;
  int len;
  double *prob;
  double (*score)[N_PEPT];
}profile;

struct PhyloTree{
  // node data
  char *label;
  double bootstrap;
  int index;
  double ewt;
  struct PhyloTree *parent;
  profile *prof;
  // children data
  int numChild;
  int DescLeaves;
  struct PhyloTree **child;
  double W;
  double w;
  double V;
  double v;
};

typedef struct PhyloTree PhyloTree;

void deletePhyloTree( PhyloTree *tree );

void deleteNode( PhyloTree *tree );

PhyloTree *copyPhyloTree( PhyloTree *tree );

PhyloTree *readPhyloTree( char *fileName );

PhyloTree *stringToPhyloTree( char *string );

int isBinaryPhyloTree( PhyloTree *tree );

int countPhyloTreeNodes( PhyloTree *tree );

PhyloTree **getAllLeaves( PhyloTree *tree , int *num);

void get_leaves( PhyloTree *tree, PhyloTree **leaves, int *num );

PhyloTree **getAllNodes( PhyloTree *tree , int *num);

void get_nodes( PhyloTree *tree, PhyloTree **leaves, int *num );

void calcDescLeaves(PhyloTree *tree);

int getPhyloTreeIndices( PhyloTree *tree, char **labels, int N);

void makeBinary(PhyloTree *node);

#endif
