// -*- coding: utf-8 -*-
#ifndef FEDATA_H
#define FEDATA_H

/* HVIS: du er utilfreds med lowercase, så husk at emacs M-l giver lovercase.
   Af interesse har også M-u, M-c */

/* http://www.cplusplus.com/doc/tutorial/classes/ */
class FEMclass{
private:
  void tot_nodes();
  void set_gll(double*, double*, int);
public:
  /* Number of element in axis direction */
  int nelx,nely,nelz;
  /* Total number of elements */
  int ne;
  /* Element lenghts in the axis directions */
  double lx,ly,lz;
  /* Number of GLL nodes per element in each axis direction */
  int ngll;
  /* Total number of nodes in directions */
  int nnx,nny,nnz;
  int nn;
  /* Constructor.  Set ne, length, gl points */
  FEMclass (int,int,int,double,double,double,int);

  /* Element numbering according to convention (ZERO INDEXED) */
  /* (i,j,k) -> (i+j*NELX+k*NELX*NELY) */
  int e(int i,int j,int k) {return (i+j*nelx+k*nelx*nely);}

  /* Node numbering according to convention (ZERO INDEXED) */
  /* (i,j,k) -> (i+j*NNX+k*NNX*NNY) */
  int n(int i,int j,int k) {return (i+j*nnx+k*nnx*nny);}

  /* Get first node in element */
  int N0(int elx,int ely, int elz) {
    return elx*(ngll-1) + ely*(ngll-1)*nnx + elz*(ngll-1)*nnx*nny;}

  /* gll points and weight */
  double *gll, *w;
};

class MATPROPclass{
public:
  /* Material properties: */
  double e, nu, thk, rho;
  MATPROPclass (double,double,double,double);
  double mu, vs;
};


/* typedef struct{ */
/*   /\* Material properties: *\/ */
/*   double e, nu, mu, thk, rho; */
/* }MATPROPtype; */
// extern MATPROPtype mat;



extern FEMclass mesh;
extern MATPROPclass mat;


#endif
