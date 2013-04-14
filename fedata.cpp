// -*- coding: utf-8 -*-

#include <stdio.h>
#include "fedata.h"
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;

/* constructor */
FEMclass::FEMclass(int i1, int i2, int i3,
                   double d1, double d2, double d3,
                   double d4, double d5,
                   int i4, int i5){

  nelx = i1, nely = i2, nelz = i3;
  lx = d1, ly = d2, lz = d3;
  offsetX = d4, offsetY = d5;
  ngll = i4; nInt = i5;
  /* set other variables */
  tot_nodes();
  set_gll();
  set_coor();
  set_Ke();
  set_edof();

  /* no dofs yet */
  numFixedDofs = 0;
  fixedDofs = new int[nn*3];
  memset(fixedDofs,0,sizeof(int)*nn*3);
}

void FEMclass::set_Ke(){
  /* must be called after tot_nodes => nen is declared */
  sizeB = new int[2];
  sizeB[0] = 6;
  sizeB[1] = nen*3;
  sizeC = new int[2];
  sizeC[0] = 6;
  sizeC[1] = 6;
}

void FEMclass::tot_nodes(){
  /* This is not really needed. Could just as well be copied to the
     constructor. But... lets keep it for the fun of it. */
  ne = nelx*nely*nelz;
  nnx = nelx*ngll-(nelx-1);
  nny = nely*ngll-(nely-1);
  nnz = nelz*ngll-(nelz-1);
  nn = nnx*nny*nnz;
  nen = pow(ngll,3);
}


void FEMclass::set_coor(){
  /* Calculates x,y,z coor from gll points */

  x=new double[nnx]; y=new double[nny]; z=new double[nnz];
  /* length of element */
  dlx = lx/nelx;
  dly = ly/nely;
  dlz = lz/nelz;

  cout << "offsetX: " << offsetX << ", offsetY: " << offsetY << endl;
  //printf("size x: %i \n",nnx);
  set_coor_helper(x,dlx,nelx,offsetX);
  set_coor_helper(y,dly,nely,offsetY);
  set_coor_helper(z,dlz,nelz,0.0);
}

void FEMclass::set_coor_helper(double *xx, double dxx,int nelxx, double offset){

  int idx=0;
  /* Helt ude at skide. Kører igennem al for mange elementer */
  for(int i=0;i<nelxx;i++){
    for(int j=0;j<ngll-1;j++){
      /* X(e+1)+X(e) = dx*(2*i+1), where X(e) is the global start coordinate of
         element e */
      // printf("idx: %i \t",idx);
      xx[idx] = 0.5*dxx*(2*i+1) + 0.5*dxx*gll[j] + offset;
      idx += 1;
    }
  }
  /* add last nodes coor */
  xx[idx] = dxx*nelxx + offset;// = lxx
}

void FEMclass::set_gll(){// (double *gll,double *w,int ngll){
  /* get gll-points and weights*/

  /* gll points for plotting */
  gll=new double[ngll];
  double *tmp = new double[ngll];
  /* numerical integration points and weights */
  intPoints = new double[nInt];
  w=new double[nInt];

  get_points(ngll,gll,tmp);
  get_points(nInt,intPoints,w);

  delete[] tmp;
}


void FEMclass::get_points(int npoints,double *point, double *weight){
/* POINT points and weights are symmetric around 0, thus only the positive are
   given explicitly in the table. The negative values are set below */
switch(npoints){
 case(2): /* 2 point : { -1 1 } (trapezoidal rule) */
   point[1] = 1;
   weight[1] = 1;
   break;
 case 3:  /* 3 point : { -1  0  1 }, weights { 1/3 4/3 1/3 } (Simpson's rule) */
   point[1] = 0.L; weight[1] =1.333333333333333333333333333333L;
   break;
 case 4:
   point[2] = 4.472135954999580e-01L; weight[2] = 8.33333333333333333333333333333e-01L;
   break;
 case 5:
   point[2] = 0.L;                    weight[2] = 7.11111111111111111111111111111e-01L;
   point[3] = 6.546536707079772e-01L; weight[3] = 5.44444444444444444444444444444e-01L;
   break;
 case 6:
   point[3] = 2.852315164806451e-01L; weight[3] = 5.548583770354862e-01L;
   point[4] = 7.650553239294655e-01L; weight[4] = 3.784749562978474e-01L;
   break;
 case 7:
   point[3] = 0.L;                    weight[3] = 4.876190476190476e-01L;
   point[4] = 4.688487934707141e-01L; weight[4] = 4.317453812098623e-01L;
   point[5] = 8.302238962785669e-01L; weight[5] = 2.768260473615680e-01L;
   break;
 case 8:
   point[4] = 2.092992179024789e-01L; weight[4] = 4.124587946587041e-01L;
   point[5] = 5.917001814331432e-01L; weight[5] = 3.411226924835035e-01L;
   point[6] = 8.717401485096070e-01L; weight[6] = 2.107042271435098e-01L;
   break;
 }
/* From
   http://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Lobatto_rules */
if (npoints>2)
  weight[npoints-1] = 2.0/(npoints*(npoints-1.0));
point[npoints-1] = 1;

/* set the negative values(weights are always positive) */
/* Notice that when npoints is odd, the center point is zero. This value should
   not be transferred */
int n = npoints/2;
if(npoints % 2  == 1){
  /* x is odd */
  n = (npoints-1)/2;
 }
for(int i=0;i<n;i++){
  point[i] = -point[npoints-1-i];
  weight[i] = w[npoints-1-i];
 }
}

/* SELECT CASE ( ng ) */
/*  case (1) */
/*   w = 2.0 */
/*   xi = 0.0 */
/*   eta =xi */
/*  case (2) */
/*   w = 1.0 */
/*   xi(1) = -3.0**(-0.5) */
/*   xi(2) = 3.0**(-0.5) */
/*   eta = xi */
/*  case (3) */
/*   w(1) = 5.0/9.0 */
/*   w(2) = 8.0/9.0 */
/*   w(3) = 5.0/9.0 */
/*   xi(1) = -0.6**(0.5) */
/*   xi(2) = 0.0 */
/*   xi(3) = 0.6**(0.5) */
/*   eta = xi */
/*  case (4)   */
/*   w(1) = (18.0+dsqrt(30d0))/36.0 */
/*   w(2) = (18.0+dsqrt(30d0))/36.0 */
/*   w(3) = (18.0-dsqrt(30d0))/36.0 */
/*   w(4) = (18.0-dsqrt(30d0))/36.0 */
/*   xi(1) = dsqrt((3d0-2d0*dsqrt(6d0/5d0))/7d0) */
/*   xi(2) = -dsqrt((3d0-2d0*dsqrt(6d0/5d0))/7d0) */
/*   xi(3) = dsqrt((3d0+2d0*dsqrt(6d0/5d0))/7d0) */
/*   xi(4) = -dsqrt((3d0+2d0*dsqrt(6d0/5d0))/7d0) */
/*   eta = xi */
// http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/#gauss_quadrature_abscissas_table
void FEMclass::set_edof(){

  /* initialize edof mat */
  edof = new int[nen*3*ne];

  int ii = 0; // element counter
  int jj = 0; // node couter
  int node0; // global number of element "0-node"
  int node; // considered element node

  /* loop elements */
  for(int elx=0; elx<nelx; elx++){
    for(int ely=0; ely<nely; ely++){
      for(int elz=0; elz<nelz; elz++){

        /* get element "0-node" */
        node0 = N0(elx,ely,elz);

        /* loop element nodes */
        jj = 0;
        for(int k=0; k<ngll; k++){
          for(int j=0; j<ngll; j++){
            for(int i=0; i<ngll; i++){

              node = node0+i+j*nnx+k*nnx*nny; // global node number

              edof[(jj+0)*ne+ii] = node*3+0; // global x-dof number
              edof[(jj+1)*ne+ii] = node*3+1; // global y-dof number
              edof[(jj+2)*ne+ii] = node*3+2; // global z-dof nubmer

              jj = jj+3;

            }
          }
        } /* end node loop */

        ii=ii+1;

      }
    }
  } /* end element loops */

}



MATPROPclass::MATPROPclass(double d1,double d2,double d3,double d4){
  /* Initialize material properties */
  e = d1, nu = d2, thk = d3, rho = d4;
  mu = e/(2*(1+nu));
  vs = sqrt(mu/rho);
}
