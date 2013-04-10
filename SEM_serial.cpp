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

  FEMclass mesh(10,10,2,10,10,2,3,4);
  /* e, nu, thk, rho; */
  MATPROPclass mat(10000,0.0,1,1);

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
  double CLF = 0.6;
  
  double dt = CLF*mesh.dlx/(mesh.ngll+1)/mat.vs;
  dt *= 0.5;
  cout << "dt: " << dt << endl;
  // CFL = 0.6; % stability number = CFL_1D / sqrt(2);
  // dt = CFL*min([lx,ly,lz])/(NGLL+1)/vs;
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

  /* cout << "M-vec: " << endl; */
  /* for(int ii=0;ii<10;ii++){ */

  /*     cout << Me[ii] << "\t"; */
  /* } */
  /*   cout << endl; */

  int tmp = sizeKe;
  cout << "K-mat: " << endl;
  for(int ii=0;ii<10;ii++){
    for(int jj=0;jj<8;jj++){
      cout << Ke[ii*tmp+jj] << "\t";
    }
    cout << endl;
  }
  //  return 1;

  /* just for testing */
  /* for(int i=0;i<1;i++) */
  /*   print_vtk(mesh,i+1,d); */

  /**********************************/
  /* STEP 3: Explicit time stepping */
  /**********************************/
  int NT = 1000;

  /* Create times step object */

  TransientSolver ExplicitSolver(&mesh,dt,NT);
  ExplicitSolver.Solve(ForceUpdater,Ke,Ce,Me);




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

