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
  e = d1, nu = d2, thk = d3, rho = d4;
  mu = e/(2*(1+nu));
  vs = sqrt(mu/rho);
}
