// -*- coding: utf-8 -*-
#ifndef FEDATA_H
#define FEDATA_H

/* HVIS: du er utilfreds med lowercase, så husk at emacs M-l giver lovercase.
   Af interesse har også M-u, M-c */

/* http://www.cplusplus.com/doc/tutorial/classes/ */
class FEMclass{
private:
  void tot_nodes();
  /* Hør! Det er da vidst ikke nødvendigt at give interne variable ind. Spørg
     Henrik */
  void set_gll();
  void get_points(int npoints,double *gll, double *w);
  void set_intPoints();
  void set_coor();
  void set_coor_helper(double *x, double dx,int nelx, double offset);
  void set_Ke();

  double offsetX, offsetY;
public:
  /* Number of element in axis direction */
  int nelx,nely,nelz;
  /* Total number of elements */
  int ne;
  /* Structure lengths in the axis directions */
  double lx,ly,lz;
  /* Element lengths in the axis direction */
  double dlx,dly,dlz;
  /* Number of GLL nodes per element in each axis direction */
  int ngll;
  /* number of nodes in one element */
  int nen;
  /* Total number of nodes in directions */
  int nnx,nny,nnz;
  int nn;

  /* Constructor.  Set ne, length, offset, gl points, int points */
  FEMclass (int,int,int,double,double,double,double,double,int,int);

  /* Element numbering according to convention (ZERO INDEXED) */
  /* (i,j,k) -> (i+j*NELX+k*NELX*NELY) */
  int e(int i,int j,int k) {return (i+j*nelx+k*nelx*nely);}

  /* Node numbering according to convention (ZERO INDEXED) */
  /* (i,j,k) -> (i+j*NNX+k*NNX*NNY) */
  int n(int i,int j,int k) {return (i+j*nnx+k*nnx*nny);}

  /* Get first node in element */
  int N0(int elx,int ely, int elz) {
    return elx*(ngll-1) + ely*(ngll-1)*nnx + elz*(ngll-1)*nnx*nny;}

  /* gll points for plotting(eq. get coordinates) */
  double *gll;
  /* Numerical integration points */
  double *intPoints, * w;
  int nInt;

  /* coordinates */
  double *x, *y, *z;

  /* array of edof indices: edof[d*ne+j], d = dof [0;nen*3], e = element [0,ne] */
  int *edof;
  void set_edof();

  /* fixed dofs */
  int *fixedDofs;
  int numFixedDofs;


  /* size of strain/displacement matrix(B) */
  int *sizeB;
  //double *B;
  /* size of 'elastic constitutive matrix'(C)
     http://en.wikipedia.org/wiki/Hooke's_law#Isotropic_materials*/
  int *sizeC;
  //double *C;
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
