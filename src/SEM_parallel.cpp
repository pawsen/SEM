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

int main(int argc, char *argv[]) {

  /* Global dimensions of mesh */
  int gnelx=47,gnely=47,gnelz=2;
  double glx=10,gly=10,glz=2;


  /***********************************/
  /* STEP 0: Initialize and setup MPI */
  /***********************************/

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
  int ngll=4, gaussPoints = 4;
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
      M[dof] += Me[n];
    }
  }
  delete[] Me;

#ifdef MY_MPI
  /* Transfer mass to neighboring nodes */
  mpi.communicate(M,false);
#endif

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

  /* Counter, used to see if each process is doing the same number of operations */
  unsigned long int count = 0;

  /* Newmark parameters */
  double gamma = 0.5;

  /* Determine time step size */
  double CLF  = 0.6;
  double dt   = CLF*mesh.dlx/(mesh.ngll+1)/mat.vs;
  dt         *= 0.5;
  double T    = 2*Ft0; /* Total integration time */
  int NT      = ceil(T/dt);

  /* Initialize kinematic fields and force vector... */
  double *d      = new double[(mesh.nn)*3];
  double *dtilde = new double[(mesh.nn)*3];
  double *v      = new double[(mesh.nn)*3];
  double *vtilde = new double[(mesh.nn)*3];
  double *a      = new double[(mesh.nn)*3];
  double *f      = new double[(mesh.nn)*3];

  /* ... and set to 0 */
  memset(d,0,sizeof(double)*mesh.nn*3);
  memset(dtilde,0,sizeof(double)*mesh.nn*3);
  memset(v,0,sizeof(double)*mesh.nn*3);
  memset(vtilde,0,sizeof(double)*mesh.nn*3);
  memset(a,0,sizeof(double)*mesh.nn*3);
  memset(f,0,sizeof(double)*mesh.nn*3);

  /* Determine center node */
  int nx_half = floor(mesh.nnx/2);
  int ny_half = floor(mesh.nny/2);
  int nz_half = floor(mesh.nnz/2);
  int center_node = mesh.n(nx_half,ny_half,nz_half);
  /* How often to save the plot? */
  int nplot = -1; double t_plot = 0.05;
  bool save_plot = true;

#ifdef MY_MPI
  /* Overlapping communications */
  MPI_Request recReqs[8], sendReqs[8];

  /* MPI-timing */
  MPI_Barrier(MPI_COMM_WORLD);
  double t1,t2;
  if(mpi.rank == 0){
    t1 = MPI_Wtime();
    printf("Number of steps, NT: %d \n",NT);
    printf("Number of real steps, NT: %d \n",(int)ceil(T/dt));
  }
#endif


  /************************/
  /* Time stepping begins */
  /************************/
  double t = 0;
  for(int it=0; it<NT; it++){
    t = t + dt;

    /* Update force */
#ifdef MY_MPI
    /* Find center domain */
    if( mpi.myCoords[0]==floor(mpi.dims[0]/2) && mpi.myCoords[1]==floor(mpi.dims[1]/2) ){
      /* Apply force in z-direction at center node */
      NodalForce(f,center_node,2,ricker(t,Ff0,Ft0));
    }
#else
    /* Apply force in z-direction at center node */
    NodalForce(f,center_node,2,ricker(t,Ff0,Ft0));
#endif

    /* Elementwise Explicit Newmark step */
#ifdef MY_OVERLAP
#include "overlap.cpp"
#else
    /* Also for serial execution */
#include "noOverlap.cpp"
#endif

    if((mpi.rank==0)&&(it%100==0))
      cout << "step:" << it << endl;

    if(save_plot){
#ifdef MY_MPI
    /* Save data.... */
    if(it % (int)ceil(t_plot/dt) == 0 ){nplot++; print_vtk(&mesh,&mpi,nplot,d);}
#else
    if(it % (int)ceil(t_plot/dt) == 0 ){nplot++; print_vtk_serial(&mesh,nplot,d);}
#endif
    } /* end save */
  }  /* end time loop */

#ifdef MY_MPI
  MPI_Barrier(MPI_COMM_WORLD);

  /* Register final time */
  if(mpi.rank==0){
    /* write .pvd-file: connect all time step */
    write_pvd(t_plot,nplot);
    t2 = MPI_Wtime();

  }

  /* SUM counters on rank 0 of global communicator*/
  unsigned long int countTotal=0;
  MPI_Reduce(&count,&countTotal,1,MPI::UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);

  /* Set time-filename according to jobname when running on HPC*/
  char filename1[30];
  if (argc != 2){
    /* A filename was not passed to the program */
    sprintf(filename1,"");
  }else
    strcpy(filename1, argv[1]);

  /* Write output */
  if(mpi.rank==0){
    char filename[30];
    sprintf(filename,"stout/timing%s.txt",filename1);

    FILE * pFile; // print time to file
    pFile = fopen (filename,"a"); // append
    fprintf (pFile, "%i \t %i \t %i \t %i \t %i \t %e \t %lu \n",mpi.size,gnelx,gnely,gnelz,NT,t2-t1,countTotal);
    printf ( "nProcs: %i, gnelx: %i, gnely: %i, gnelz: %i , NT: %i, Time: %e, totCount: %lu \n",mpi.size,gnelx,gnely,gnelz,NT,t2-t1,countTotal);
    fclose (pFile);

  }

  /* /\* Print the number of "operations" of this proc to file *\/ */
  /* FILE * pFile1; // print time to file */
  /* char fName[100]; */
  /* sprintf(fName,"stout/numOp_rank%i%s.txt",mpi.rank,filename1); */
  /* pFile1 = fopen (fName,"a"); // append */
  /* fprintf (pFile1, "%i \t %i \t %lu \n",mpi.rank,mpi.size,count); */
  /* printf ( "myRank: %i, nProcs: %i, count: %lu \n",mpi.rank,mpi.size,count); */
  /* fclose (pFile1); */

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

/* Used to to apply a nodal force */
void NodalForce(double* F, int node, int dof, double val){
  // *F = pointer to force array
  // node = global node number
  // dof = 0 1 2
  // val = force value
  F[ node*3+dof ] = val;

}

double min3(double x, double y, double z){
  if (x < y)
    if (x < z) return x;
  if (y < z) return y;
  return z;
}
