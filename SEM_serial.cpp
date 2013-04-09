// -*- coding: utf-8 -*-
/* FEM implementation using spectral elements */
/* Lagrange polynomials are used as shape functions.

   Since the integration points are gauss Lobatto points, and Lagrange
   polonomials are orthogonal in these points(and only in these points(η,ζ)),
   the integration over the "mass" yields(2D)
   M_{ij}=∫₁¹ ρ⋅L(η)⋅L(ζ) dη dζ = ρ⋅δ_{ij}
   eg. the mass matrix is diagonal

   http://www.cs.berkeley.edu/~demmel/cs267/
*/


/* Use extern C to tell the C++ compiler that the included library is compiled
  with C and not C++
  http://stackoverflow.com/a/67930/1121523 */
extern "C" {
#include <cblas.h>
  //#include <atlas/cblas.h>
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for memset */
#include <math.h>
#include "fedata.h"
#include "vtk.h"
#include "aux.h"
#include "transient.h"
#include <assert.h>
#include <typeinfo>

#include <iostream>
using namespace std;

/* Force function */
double ricker(double t, double f0, double to);

/* Set nodal force */
void NodalForce(double* F, int node, int dof, double val);

void ForceUpdater(FEMclass* mesh, double* F, double t);

int main(){

  /****************************************/
  /* STEP 1: Initialize mesh and material */
  /****************************************/

  FEMclass mesh(2,2,10,1,1,1,3);
  /* e, nu, thk, rho; */
  MATPROPclass mat(10000,0,1,1);

  /* cout << "nodes x: "; */
  /* for(int i=0;i<mesh.nnx;i++){ */
  /*   cout << mesh.x[i] << '\t'; */
  /* } */
  /* cout << endl; */


  /**************************/
  /* STEP 2: INITIALIZATION */
  /**************************/

  /* Force function. Gauss like: Mexican hat */
  double Ff0 = 0.25;            /* fundamental frequency */
  double Ft0 = 1.5/Ff0;
  double Ft;

  /* ############################################ */
  /* Testet output HERTIL. Giver samme som SEM.m  */

  /* The timestep dt is set by the stability condition CLF */
  /* http://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition#The_two_and_general_n-dimensional_case */

  /* Ask Boyan  CLF and Δx/Δy Is that the element length, or minimum distance
     between nodes?. */
  double CLF = 1;
  double dt = 0.0000001;
  double half_dt = 0.5*dt;

  /* Global elements first node(N0). Used to calculate rest of elements nodes*/
  double No, node;
  //int edof[sizeKe];
  int centerNode = mesh.n(floor(mesh.nnx/2),floor(mesh.nny/2),floor(mesh.nnz/2)); /* force applied here */

  /* allocate kinematic fields and forces */
  double *d = new double[mesh.nn*3];
  double *v = new double[mesh.nn*3];
  double *a = new double[mesh.nn*3];
  double *f = new double[mesh.nn*3];

  /* set entries to zero */
  memset(d,0, sizeof d), memset(v,0, sizeof v), memset(a,0, sizeof a);
  memset(f,0,sizeof f);


  /******************************/
  /* Construct element matrices */
  /******************************/

  /* size along each dimension (number of element dofs) */
  int sizeKe = mesh.nen*3;
  /* http://www.cplusplus.com/reference/cstring/memset/ */
  /* allocate dynamically. Otherwise there's stack overflow for ngll>6 */
  double *Ke = new double[sizeKe*sizeKe];
  double Me[sizeKe];
  double Ce[sizeKe];
  memset(Ke,0,sizeof(Ke));
  memset(Me,0,sizeof(Me));
  memset(Ce,0,sizeof(Ce));

  matmul_test();
  construct_Ke(&mat,&mesh, Ke);
  construct_Me(&mat,&mesh, Me);

  /* MASS PROP. DAMPING:
     C=αM+βK, ξ=0.5(α/ω+βω)
     Choosing ξ=0.01 at ω=2π*f0
     -> α=4πξ*f0*/
  double xi = 0.01;              /* 1% damping at selected frequency */
  double alpha = 4*M_PI*xi*Ff0;  /* Ff0 from force function */
  for(int i=0; i<mesh.nen*3; i++)
    Ce[i] = alpha*Me[i];

  /* just for testing */
  for(int i=0;i<100;i++)
    print_vtk(mesh,i+1,d);

  /**********************************/
  /* STEP 3: Explicit time stepping */
  /**********************************/
  int NT = 2;

  /* Create times step object */

  /**********************************************************************/
  /* SKAL DU IKKE ALLOKERE d,v,a osv? Så tror jeg det fungerer bedre :) */
  /**********************************************************************/
  /* TransientSolver ExplicitSolver(&mesh,dt,NT); */
  /* ExplicitSolver.Solve(ForceUpdater,Ke,Ce,Me); */
  

}

double ricker(double t, double f0, double t0){
  /* Return the force, given by an: */
  /* Mexican hat wavelet */
  /* http://en.wikipedia.org/wiki/Mexican_hat_wavelet */

  double arg = M_PI*f0*(t-t0);
  arg = arg*arg;
  double f = (2*arg-1)*exp(-arg);
  return f;

}

/* Apply nodal force function */
void NodalForce(double* F, int node, int dof, double val){
  // *F = pointer to force array
  // node = global node number
  // dof = 0 1 2
  // val = force value
  F[ node*3+dof ] = val;
}

/* Force function */
void ForceUpdater(FEMclass* mesh,double* F,double t){
  
  /* det er selvfølgeligt lidt uheldigt at det skal hardcodes her igen.. */
  double Ff0 = 0.25;
  double Ft0 = 1.5/Ff0;

  /* center node */
  int nx_half = floor(mesh->nnx/2);
  int ny_half = floor(mesh->nny/2);
  int nz_half = floor(mesh->nnz/2);
  int center_node = mesh->n(nx_half,ny_half,nz_half);

  /* Apply force in z-direction */
  NodalForce(F,center_node,2,ricker(t,Ff0,Ft0));

}

