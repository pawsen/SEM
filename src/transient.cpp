// -*- coding: utf-8 -*-

#include "transient.h"
#include "vtk.h"
#include <stdio.h>
#include <string.h> /* for memset */
#include <iostream>
#include <math.h>
using namespace std;

TransientSolver::TransientSolver(FEMclass *mesh_in, double dt_in=0.0, int NT_in = 0)
{
  /* default parameters */
  beta = 0.0;
  gamma = 0.5;


  mesh = mesh_in;
  dt = dt_in;
  NT = NT_in;

  /* Allocate kinematic fields and forces */
  d = new double[(mesh->nn)*3];
  dtilde = new double[(mesh->nn)*3];
  v = new double[(mesh->nn)*3];
  vtilde = new double[(mesh->nn)*3];
  a = new double[(mesh->nn)*3];
  f = new double[(mesh->nn)*3];

  /* Set to 0 */

  memset(d,0,sizeof(double)*mesh->nn*3);
  memset(dtilde,0,sizeof(double)*mesh->nn*3);
  memset(v,0,sizeof(double)*mesh->nn*3);
  memset(vtilde,0,sizeof(double)*mesh->nn*3);
  memset(a,0,sizeof(double)*mesh->nn*3);
  memset(f,0,sizeof(double)*mesh->nn*3);

  /* for(int i=0;i < mesh->nn*3; i++) */
  /*   cout << "d[0]: " << d[i] << endl; */
};

void TransientSolver::ExplicitCDSStep(const double *Ke, const double *KSpring,
                                      const double *C, const double *M){

  /*********************************************************************/
  /* Newmark Time Integration overview                                 */
  /*                                                                   */
  /* ( Described in 16.5.1 of Finite Element Method by Wiberg et al.)  */
  /*                                                                   */
  /* METHOD                TYPE      γ    β    STABILITY               */
  /* Average acceleration  implicit  1/2  1/4  unconditional           */
  /* Linear acceleration   implicit  1/2  1/6  Ωcr = 2*3^(1/2)         */
  /* Fox-Goodwin           implicit  1/2  1/12 Ωcr = 6^(1/2)           */
  /* Central Difference    explicit  1/2  0    Ωcr = 2    <--THIS ONE! */
  /*                                                                   */
  /* This algorithm assumes diagonal mass- and damping matrix             */
  /*********************************************************************/


  //cout << "a \n";



  //cout << "b \n";

  /* STEP 1: Define Predictors */
  for(int i=0; i<(mesh->nn)*3; i++){
    dtilde[i] = d[i] + v[i]*dt + 0.5*dt*dt*(1-2*beta)*a[i];
    vtilde[i] = v[i] + (1-gamma)*dt*a[i];

    a[i] = 0.0; // while looping, also reset acceleration for step 2
  }

/* "backward" communications makes it explode */
#ifdef MY_MPI
  mpi->communicate(dtilde,true);
  mpi->communicate(vtilde,true);
#endif


  /* STEP 2: Solve for Acceleration */

  /* loop elements, e */
  /* NOTE: OVERVEJ: */
  /* for(int e=0; e<mesh.ne; e++){ */
  int dof; // single entry of edof-vector
  double prod = 0; // single entry of ke*dtilde product
  int e = 0;
  for(int elx=0; elx<mesh->nelx; elx++){
    for(int ely=0; ely<mesh->nely; ely++){
      for(int elz=0; elz<mesh->nelz; elz++){

        /* loop dofs, i in edofs-vec for this element */
        for(int i=0; i<mesh->nen*3; i++){

          /* prod = dof-th entry of vector from ke*dtilde product */
          prod = 0;
          for(int j=0; j<mesh->nen*3; j++){
            prod = prod + Ke[i*mesh->nen*3+j]*dtilde[ mesh->edof[j*mesh->ne+e] ];
          }

          /* cout << "prod done, i,e=" << i<< ","<<e << "\n"; */

          // global index of the considered dof of the element
          dof = mesh->edof[i*mesh->ne+e];
          /* cout << mesh->edof[i*mesh->ne+e] << " hej " <<"\t"; */
          /* cout << dof << "\t"; */
          
          /* RHS, Division cant be done on element basis because M^-1 ≠ ∑Me^-1 */
          a[dof] = a[dof] + ( f[dof] - C[dof]*vtilde[dof] - prod - KSpring[dof]*dtilde[dof] );
	  /* a[dof] = a[dof] + ( f[dof] - C[dof]*vtilde[dof] - prod - KSpring[dof]*dtilde[dof] )/(M[i] + gamma*dt*C[i]); */

        } /* end dof loop */

        e++; /* increment element */

      }
    }
  } /* end of element loops */

#ifdef MY_MPI
  mpi->communicate(a,false); /* "forward" communication */
#endif

  /* Divide with global mass vector */
  for(int i=0;i<mesh->nn*3;i++){
    a[i] /= (M[i] + gamma*dt*C[i]);
  }

/* #ifdef MY_MPI */
/*   mpi->communicate(a,false); /\* "forward" communication *\/ */
/* #endif */

  /* STEP 3: Calculate Displacement and Velocity as Correctors */
  for(int i=0; i<mesh->nn*3; i++){
    d[i] = dtilde[i] + beta*dt*dt*a[i];
    v[i] = vtilde[i] + gamma*dt*a[i];
  }

}

void TransientSolver::Solve(void (*ft)(FEMclass*,double*,double,int),
                            const double *Ke,const double *KSpring,
                            const double *C, const double *M){
 
  /* Time loop */
  double t = 0;
  /* How often to save the plot? */
  int nplot = -1; double t_plot = 0.05;
  for(int it=0; it<NT; it++){
    t = t + dt;

    /* Update force */
#ifdef MY_MPI
    if( mpi->myCoords[0]==floor(mpi->dims[0]/2) && mpi->myCoords[1]==floor(mpi->dims[1]/2) ){
      (*ft)(mesh,f,t,it);
    }
#else
    (*ft)(mesh,f,t,it);
#endif

    /* Explicit Newmark step */
    ExplicitCDSStep(Ke,KSpring,C,M);

#ifdef MY_MPI
    /* Save data.... */
    //if(it % (int)ceil(t_plot/dt) == 0 ){nplot++; print_vtk(mesh,mpi,nplot,d);}

    /* if( mpi->rank==1 )// mpi->myCoords[0]==floor(mpi->dims[0]/2) && mpi->myCoords[1]==floor(mpi->dims[1]/2) ) */
    /*   print_vtk(mesh,it,d); */

    /* Write out info every 100 steps */
    if( (mpi->rank == 0) && (it % 100 == 0) ){
      printf("step=%d\n",it);
      cout << "rank: "<< mpi->rank << ", t: " << t << "\t" << "dt: " << dt << ", d[0]=" << d[0] << endl;
    }
#else
    /* Save data.... */
    if(it % (int)ceil(t_plot/dt) == 0 ){nplot++; print_vtk_serial(mesh,nplot,d);}

    if (it % 100 == 0)
      cout << "t: " << t << "\t" << "dt: " << dt << ", d[0]=" << d[0] << endl;
#endif

  } /* end time loop */


#ifdef MY_MPI
  /* write .pvd-file: connect all time step */
  write_pvd(t_plot,nplot);
#endif

}

#ifdef MY_MPI
void TransientSolver::transferMPI(MPIClass *mpi_in){
  /* Make the mpi-class accessible for communicating between nodes */
  mpi = mpi_in;
}
#endif
