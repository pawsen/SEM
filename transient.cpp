// -*- coding: utf-8 -*-

#include "transient.h"
#include "vtk.h"
#include <stdio.h>
#include <string.h> /* for memset */
#include <iostream>
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
  memset(d,0, sizeof d);
  memset(dtilde,0, sizeof dtilde);
  memset(v,0, sizeof v);
  memset(vtilde,0, sizeof vtilde);
  memset(a,0, sizeof a);
  memset(f,0,sizeof f);
};

void TransientSolver::ExplicitCDSStep(const double *Ke,
                                      const double *Ce, const double *Me){

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

  int dof; // single entry of edof-vector
  double prod = 0; // single entry of ke*dtilde product

  //cout << "b \n";

  /* STEP 1: Define Predictors */
  for(int i=0; i<(mesh->nn)*3; i++){
    dtilde[i] = d[i] + v[i]*dt + 0.5*dt*dt*(1-2*beta)*a[i];
    vtilde[i] = v[i] + (1-gamma)*dt*a[i];

    a[i] = 0.0; // while looping, also reset acceleration for step 2
  }

  //cout << "step1 done" << endl;

  /* STEP 2: Solve for Acceleration */

  /* cout << "dof : \n"; */

  /*  int ee = 0; */
  /* for(int elx=0; elx<mesh->nelx; elx++){ */
  /*   for(int ely=0; ely<mesh->nely; ely++){ */
  /*     for(int elz=0; elz<mesh->nelz; elz++){ */
  /*       for(int i=0; i<mesh->nen*3; i++){ */
  /*         //dof = mesh->edof[i*mesh->nen*3+ee]; */
  /*         //cout << mesh->edof[i*mesh->nen*3+ee]  <<"\t"; */
  /*      dof = mesh->edof[i*mesh->ne+ee]; */
  /*         cout << "i,e=" << i << "," << ee << "," << mesh->edof[i*mesh->ne+ee]  <<"\n"; */
  /*       } */
  /*     ee++; */
  /*    return; */
  /*     } */
  /*   } */
  /* } */
  /* cout << "\n ------- DONE --------- \n"; */

  /* loop elements, e */
  int e = 0;
  for(int elx=0; elx<mesh->nelx; elx++){
    for(int ely=0; ely<mesh->nely; ely++){
      for(int elz=0; elz<mesh->nelz; elz++){

        /* loop dofs, i in edofs-vec for this element */
        for(int i=0; i<mesh->nen*3; i++){

          /* prod = dof-th entry of vector from ke*dtilde product */
          prod = 0;
          for(int j=0; j<mesh->nen*3; j++){
            prod = prod + Ke[i*mesh->nen*3+j]*dtilde[ mesh->edof[j*mesh->ne+e]];
          }

          /* cout << "prod done, i,e=" << i<< ","<<e << "\n"; */

          // global index of the considered dof of the element
          dof = mesh->edof[i*mesh->ne+e];
          /* cout << mesh->edof[i*mesh->ne+e] << " hej " <<"\t"; */
          /* cout << dof << "\t"; */
          // no springs yet!
          a[dof] = a[dof] + ( f[dof] - Ce[i]*vtilde[dof] - prod ) /
            (Me[i]+gamma*dt*Ce[i]);

        } /* end dof loop */

        e++; /* increment element */

      }
    }
  } /* end of element loops */

  /* enforce zero aceleration at fixed dofs */
  for(int i=0; i<mesh->numFixedDofs; i++)
    a[ mesh->fixedDofs[i] ] = 0;

  //cout << "step2 done" << endl;


  /* STEP 3: Calculate Displacement and Velocity as Correctors */
  for(int i=0; i<mesh->nn*3; i++){
    d[i] = dtilde[i] + beta*dt*dt*a[i];
    v[i] = vtilde[i] + gamma*dt*a[i];
  }

}

void TransientSolver::Solve(void (*ft)(FEMclass*,double*,double),
                            const double *Ke,const double *Ce, const double *Me){

  /* Time loop */
  int t = 0;
  for(int it=0; it<NT; it++){
    t = t + dt;

    /* Update force */
    (*ft)(mesh,f,t);

    //cout << "går ind i explicitCDSStep!" << endl;

    /* Elementwise Explicit Newmark step */
    ExplicitCDSStep(Ke,Ce,Me);

    /* Save data.... */
    print_vtk(mesh,it,d);

    /* Write out info every 100 steps */
    if(it % 1000 == 0){
      printf("step=%d\n",it);
    }

  }
  /* end time loop */

}
