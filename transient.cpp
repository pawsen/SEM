// -*- coding: utf-8 -*-

#include "transient.h"
#include <stdio.h>
#include <string.h> /* for memset */
TransientSolver::TransientSolver(FEMclass *mesh_in, double dt_in=0.0, int NT_in = 0)
{
  /* default parameters */
  beta = 0.0;
  gamma = 0.5;

  mesh = mesh_in;
  dt = dt_in;
  NT = NT_in;
};

void TransientSolver::ExplicitCDSStep(const double *F,const double *Ke,
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

  int dof; // single entry of edof-vector
  double prod = 0; // single entry of ke*dtilde product

  /* STEP 1: Define Predictors */
  for(int i=0; i<mesh->nn*3; i++){
    dtilde[i] = d[i] + v[i]*dt + 0.5*dt*dt*(1-2*beta)*a[i];
    vtilde[i] = v[i] + (1-gamma)*dt*a[i];

    a[i] = 0.0; // while looping, also reset acceleration for step 2
  }

  /* STEP 2: Solve for Acceleration */

  /* loop elements, e */
  int e = 0;
  for(int elx=0; elx<mesh->nelx; elx++){
    for(int ely=0; ely<mesh->nely; ely++){
      for(int elz=0; elz<mesh->nelz; elz++){
        e++;

        /* loop dofs, i in edofs-vec for this element */
        for(int i=0; i<mesh->nen*3; i++){

          /* prod = dof-th entry of vector from ke*dtilde product */
          prod = 0;
          for(int j=0; j<mesh->nen*3; j++){
            prod = prod + Ke[i*mesh->nen*3+j]*dtilde[ mesh->edof[j*mesh->nen*3+e]  ];
          }

          // global index of the considered dof of the element
          dof = mesh->edof[i*mesh->nen*3+e];
          // no springs yet!
          a[dof] = a[dof] + ( F[dof] - Ce[dof]*vtilde[dof] - prod ) / (Me[dof]+gamma*dt*Ce[dof]);

        } /* end dof loop */

      }
    }
  } /* end of element loops */

  /* enforce zero aceleration at fixed dofs */
  for(int i=0; i<mesh->numFixedDofs; i++)
    a[ mesh->fixedDofs[i] ] = 0;


  /* STEP 3: Calculate Displacement and Velocity as Correctors */
  for(int i=0; i<mesh->nn*3; i++){
    d[i] = dtilde[i] + beta*dt*dt*a[i];
    v[i] = vtilde[i] + gamma*dt*a[i];
  }

}

void TransientSolver::Solve(void (*ft)(FEMclass*,double*,double),
                            const double *Ke,const double *Ce, const double *Me){

  double *F = new double[mesh->nn*3];
  memset(F,0,sizeof(F)); /* Set to 0 */

  /* Time loop */
  int t = 0;
  for(int it=0; it<NT; it++){
    t = t + dt;

    /* Update force */
    (*ft)(mesh,F,t);

    /* Elementwise Explicit Newmark step */
    ExplicitCDSStep(F,Ke,Ce,Me);

    /* Save data.... */


    /* Write out info every 100 steps */
    if(it % 1000 == 0){
      printf("step=%d",it);
    }

  }
  /* end time loop */

  delete[] F;

}
