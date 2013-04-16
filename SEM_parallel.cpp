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

/* Run as: */
/* make && mpirun -np 8 -mca btl tcp,self SEM_p */


/* Use extern C to tell the C++ compiler that the included library is compiled
  with C and not C++
  http://stackoverflow.com/a/67930/1121523 */
extern "C" {
#include <cblas.h>
  //#include <atlas/cblas.h>
}

// #define MY_MPI

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for memset */
#include <math.h>
#include "fedata.h"
#include "vtk.h"
#include "aux.h"
#include "transient.h"
#ifdef MY_MPI
#include "mpidata.h"
#include <mpi.h>
#endif

#include <iostream>
using namespace std;

double min3(double x,double y,double z);
/* Force function */
double ricker(double t, double f0, double to);

/* Set nodal force */
void NodalForce(double* F, int node, int dof, double val);
void ForceUpdater(FEMclass* mesh, double* F, double t, int it);

int main(int argc, char *argv[]) {

  /* Global dimensions of mesh */
  int gnelx=12,gnely=12,gnelz=2;
  double glx=12,gly=12,glz=2;

  /***********************************/
  /* STEP 0: Initialize and setup MPI */
  /***********************************/


  /* Dette er ikke så godt. Jeg forestiller mig at mpi klassen får  gnel* og
     gl* ind. Men ikke som herunder. I stedet via fedata. Dvs gnel* og gl* er
     på een eller anden måde initialiseret i fedata.h. Men hvordan? For den
     rigtige constructor i fedata.h kan først kaldes når mpi klassen er kørt og
     derved har bestemt lokale data.*/

  double lx,ly,lz, offsetX = 0, offsetY = 0;
  int nelx,nely,nelz;
#ifdef MY_MPI
  DataStruct data;
  data.inInt[0] = gnelx;data.inInt[1] = gnely;data.inInt[2] = gnelz;
  data.inDouble[0] = glx; data.inDouble[1] = gly; data.inDouble[2] = glz;

  /* initialize */
  MPIClass mpi(argc,argv,&data);

  lx = data.outDouble[0]; ly = data.outDouble[1];  lz = data.outDouble[2];
  nelx = data.outInt[0];  nely = data.outInt[1];  nelz = data.outInt[2];
  offsetX = mpi.offsetX; offsetY = mpi.offsetY;
#else
  lx = glx; ly = gly; lz = glz;
  nelx = gnelx; nely = gnely; nelz = gnelz;
#endif


  /**********************************************/
  /* STEP 1: Initialize local mesh and material */
  /**********************************************/
  int ngll= 5;
  int gaussPoints = 4;
  cout << "ngll: " << ngll << ", gauss points: "  << gaussPoints << endl ;
  FEMclass mesh(nelx,nely,nelz,lx,ly,lz,offsetX,offsetY,ngll,gaussPoints);
  /* e, nu, thk, rho; */
  MATPROPclass mat(10000,0.0,1,1);


#ifdef MY_MPI
  /* share mesh with mpi and init send/recv buffers*/
  mpi.initBuffers(&mesh);
#endif



  /**************************/
  /* STEP 2: INITIALIZATION */
  /**************************/

  /* Force function. Gauss like: Mexican hat */
  double Ff0 = 0.25;            /* fundamental frequency */
  double Ft0 = 1.5/Ff0;
  double Ft;


  /******************************/
  /* Construct element matrices */
  /******************************/

  /* size along each dimension (number of element dofs) */
  int sizeKe = mesh.nen*3;

  /* http://www.cplusplus.com/reference/cstring/memset/ */
  /* allocate dynamically. Otherwise there's stack overflow for ngll>6 */
  double *Ke = new double[sizeKe*sizeKe];
  double *Me = new double[sizeKe];
  memset(Ke,0,sizeof(double)*sizeKe*sizeKe);
  memset(Me,0,sizeof(double)*sizeKe);

  //matmul_test();
  construct_Ke(&mat,&mesh, Ke);
  construct_Me(&mat,&mesh, Me);

/* int tmp = sizeKe; */
/*   cout << "K-mat: " << endl; */
/*   for(int ii=0;ii<10;ii++){ */
/*     for(int jj=0;jj<8;jj++){ */
/*       cout << Ke[ii*tmp+jj] << "\t"; */
/*     } */
/*     cout << endl; */
/*   } */
/*   return 1; */

  /* Local(but for all elements on the proc) mass vector with contribution from
     neighboring nodes */
  /* Because M^-1 ≠ ∑Me^-1 */
  double *M = new double[mesh.nn*3];
  double *C = new double[mesh.nn*3];
  memset(M,0, sizeof(double)*mesh.nn*3);
  memset(C,0, sizeof(double)*mesh.nn*3);

  int dof;
  for(int e=0; e<mesh.ne; e++){
    /* loop dofs, i in edofs-vec for this element */
    for(int n=0; n<mesh.nen*3; n++){
      // global index of the considered dof of the element
      dof = mesh.edof[n*mesh.ne+e];
      /* cout << dof << "\t"; */
      M[dof] += Me[n];
    }
  }
  delete[] Me;

#ifdef MY_MPI
  /* Transfer mass to neighboring nodes */
  mpi.communicate(M,false);
#endif

  //  MPI_Finalize();
  //return 1;

  /* MASS PROP. DAMPING:
     C=αM+βK, ξ=0.5(α/ω+βω)
     Choosing ξ=0.01 at ω=2π*f0
     -> α=4πξ*f0*/
  double xi = 0.01;              /* 1% damping at selected frequency */
  double alpha = 4*M_PI*xi*Ff0;  /* Ff0 from force function */
  for(int i=0; i<mesh.nn*3; i++)
    C[i] = alpha*M[i];


  /* Add springs to all dofs at (x,y,z=0) */
  double k_spring = mat.e/100; /* default spring stiffness */
  double *KSpring = new double[mesh.nn*3]; /* global vector of spring stiffnesses */
  memset(KSpring,0,sizeof(double)*mesh.nn*3);       /* init to 0 */
  for(int i=0; i<mesh.nnx*mesh.nny*3; i++)
    KSpring[i] = k_spring;


  /**********************************/
  /* STEP 3: Explicit time stepping */
  /**********************************/
  double CLF = 0.6;

  double min_dx = min3(mesh.dlx,mesh.dly,mesh.dlz);
  double dt = CLF*min_dx/(mesh.ngll+1)/mat.vs;
  dt *= 0.5;
  /* Total integration time */
  double T = 2*Ft0;
  int NT = ceil(T/dt);
  
  // CFL = 0.6; % stability number = CFL_1D / sqrt(2);
  // dt = CFL*min([lx,ly,lz])/(NGLL+1)/vs;

  /* Create timestep object */
  TransientSolver ExplicitSolver(&mesh,dt,NT);
#ifdef MY_MPI
ExplicitSolver.transferMPI(&mpi);
#endif
 ExplicitSolver.Solve(ForceUpdater,Ke,KSpring,C,M);

#ifdef MY_MPI
  MPI_Finalize();
#endif

  return 1;

}

double ricker(double t, double f0, double t0){
  /* Return the force, given by an: */
  /* Mexican hat wavelet */
  /* http://en.wikipedia.org/wiki/Mexican_hat_wavelet */

  double arg = M_PI*f0*(t-t0);
  arg = arg*arg;
  double f = (2*arg-1.0)*exp(-arg);
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
void ForceUpdater(FEMclass* mesh,double* F,double t,int it){

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

double min3(double x, double y, double z){
  if (x < y)
    if (x < z) return x;
  if (y < z) return y;
  return z;
}
