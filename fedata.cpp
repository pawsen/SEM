// -*- coding: utf-8 -*-

#include "fedata.h"
#include <math.h>

/* constructor */
FEMclass::FEMclass(int i1, int i2, int i3,
                   double d1, double d2, double d3,
                   int i4){

  nelx = i1, nely = i2, nelz = i3;
  lx = d1, ly = d2, lz = d3;
  ngll = i4;
  /* set other variables */
  tot_nodes();
  set_gll();
}

void FEMclass::tot_nodes(){
  /* This is not really needed. Could just as well be copied to the
     constructor. But... lets keep it for the fun of it. */
  ne = nelx*nely*nelz;
  nnx = nelx*ngll-(nelx-1);
  nny = nely*ngll-(nely-1);
  nnz = nelz*ngll-(nelz-1);
  nn = nnx*nny*nnz;
}



MATPROPclass::MATPROPclass(double d1,double d2,double d3,double d4){
  /* Initialize material properties */
  e = d1, nu = d2, thk = d3, rho = d4;
  mu = e/(2*(1+nu));
  vs = sqrt(mu/rho);
}

void FEMclass::set_gll(double *gll,double *w,int ngll){
  /* get gll-points and weights*/

  gll=new double[ngll];
  w=new double[ngll];

  /* GLL points and weights are symmetric around 0, thus only the positive are
     given explicitly in the table. The negative values are set below */
  switch(ngll){
  case(2): /* 2 point : { -1 1 } (trapezoidal rule) */
    gll[1] = 1;
    w[1] = 1;
    break;
  case 3:  /* 3 point : { -1  0  1 }, weights { 1/3 4/3 1/3 } (Simpson's rule) */
    gll[1] = 0.L; w[1] =1.333333333333333333333333333333L;
  break;
  case 4:
    gll[2] = 4.472135954999580e-01L; w[2] = 8.33333333333333333333333333333e-01L;
    break;
  case 5:
    gll[2] = 0.L;                    w[2] = 7.11111111111111111111111111111e-01L;
    gll[3] = 6.546536707079772e-01L; w[3] = 5.44444444444444444444444444444e-01L;
    break;
  case 6:
    gll[3] = 2.852315164806451e-01L; w[3] = 5.548583770354862e-01L;
    gll[4] = 7.650553239294655e-01L; w[4] = 3.784749562978474e-01L;
    break;
  case 7:
    gll[3] = 0.L;                    w[3] = 4.876190476190476e-01L;
    gll[4] = 4.688487934707141e-01L; w[4] = 4.317453812098623e-01L;
    gll[5] = 8.302238962785669e-01L; w[5] = 2.768260473615680e-01L;
    break;
  case 8:
    gll[4] = 2.092992179024789e-01L; w[4] = 4.124587946587041e-01L;
    gll[5] = 5.917001814331432e-01L; w[5] = 3.411226924835035e-01L;
    gll[6] = 8.717401485096070e-01L; w[6] = 2.107042271435098e-01L;
    break;
  }
  /* From
     http://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Lobatto_rules */
  if (ngll>2)
    w[ngll-1] = 2.0/(ngll*(ngll-1.0));
  gll[ngll-1] = 1;

  /* set the negative values */
  /* Notice that when ngll is odd, the center point is zero. This value should
     not be transferred */
  int n = ngll/2;
  if(ngll % 2  == 1){
    /* x is odd */
    n = (ngll-1)/2;
  }
  for(int i=0;i<n;i++){
    gll[i] = -gll[ngll-1-i];
    w[i] = -w[ngll-1-i];
  }
}
