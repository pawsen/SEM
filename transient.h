// -*- coding: utf-8 -*-
#ifndef TRANSIENT_H
#define TRANSIENT_H

#include "fedata.h"
             
class TransientSolver{
private:

  /* solver parameters */
  double beta;
  double gamma; 
  
  /* FEM mesh */
  FEMclass *mesh;

  /* kinematic fields (predictors) */
  double *dtilde;
  double *vtilde;
  
  /* take a single time step */
  void ExplicitCDSStep(const double *Ke,const double *Ce, const double *Me);

public:

  /* kinematic fields */
  double *d;
  double *v;
  double *a;

  /* source vector */
  double *f;

  /* stepsize */
  double dt;

  /* number of time steps */
  double NT;  

  /* constructor */
  TransientSolver(FEMclass *mesh_in,double dt_in, int NT_in);

  /* time integrate t = [0-Tmax], note the first argument is a function pointer. */
  void Solve(void(*ft)(FEMclass*,double*,double),const double *Ke,const double *Ce, const double *Me);
};

#endif
