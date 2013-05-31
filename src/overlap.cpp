/* STEP 1: Define Predictors */

/* loop dofs by looping nodes */
for(int i=0; i<(mesh.nnx); i++){
  for(int j=0; j<(mesh.nny); j++){
    for(int k=0; k<(mesh.nnz); k++){

      /* d1 = dof+0, d2 = dof+1, d3 = dof+2 */
      dof = mesh.n(i,j,k)*3;

      dtilde[dof+0] = d[dof+0] + v[dof+0]*dt + 0.5*dt*dt*a[dof+0];
      dtilde[dof+1] = d[dof+1] + v[dof+1]*dt + 0.5*dt*dt*a[dof+1];
      dtilde[dof+2] = d[dof+2] + v[dof+2]*dt + 0.5*dt*dt*a[dof+2];

      vtilde[dof+0] = v[dof+0] + (1-gamma)*dt*a[dof+0];
      vtilde[dof+1] = v[dof+1] + (1-gamma)*dt*a[dof+1];
      vtilde[dof+2] = v[dof+2] + (1-gamma)*dt*a[dof+2];

      /* while looping, also reset acceleration for step 2 */
      a[dof+0] = 0.0;
      a[dof+1] = 0.0;
      a[dof+2] = 0.0;
    }
  }
 }

/* "backward" communication */
#ifdef MY_MPI
mpi.ReqComm(dtilde,true,recReqs,sendReqs,vtilde); /* overlapping "backward" communication, note that we send both vectors at the same time! */
#endif


/* STEP 2: Solve for Acceleration */

int dof; // single entry of edof-vector
double prod = 0; // single entry of ke*dtilde product

/* loop INTERNAL elements */
for(int i=1; i<(mesh.nelx-1); i++){
  for(int j=1; j<(mesh.nely-1); j++){
    for(int k=0; k<(mesh.nelz); k++){
      int e = mesh.e(i,j,k);
      /* loop dofs, i in edofs-vec for this element */
      for(int ii=0; ii<mesh.nen*3; ii++){

        /* prod = dof-th entry of vector from ke*dtilde product */
        prod = 0;
        for(int jj=0; jj<mesh.nen*3; jj++){
          count++;
          prod = prod + Ke[ii*mesh.nen*3+jj]*dtilde[ mesh.edof[jj*mesh.ne+e] ];
        }
        // global index of the considered dof of the element
        dof = mesh.edof[ii*mesh.ne+e];

        /* RHS, Division cant be done on element basis because M^-1 ≠ ∑Me^-1 */
        a[dof] = a[dof] + ( f[dof] - C[dof]*vtilde[dof] - prod -
                            KSpring[dof]*dtilde[dof] );

      } /* end dof loop */
    }
  }
 } /* end of element loops */

/* Complete receive communications */
MPI_Waitall(8,recReqs,MPI_STATUS_IGNORE);
mpi.EmptyBuffer(dtilde,true,vtilde);

/* loop BOUNDARY elements */
for(int i=0; i<(mesh.nelx); i+=(mesh.nelx)){
  for(int j=0; j<(mesh.nely); j+=(mesh.nely)){
    for(int k=0; k<(mesh.nelz); k++){

      int e = mesh.e(i,j,k);

      /* loop dofs, i in edofs-vec for this element */
      for(int ii=0; ii<mesh.nen*3; ii++){

        /* prod = dof-th entry of vector from ke*dtilde product */
        prod = 0;
        for(int jj=0; jj<mesh.nen*3; jj++){

          count++;
          prod = prod + Ke[ii*mesh.nen*3+jj]*dtilde[ mesh.edof[jj*mesh.ne+e] ];
        }

        // global index of the considered dof of the element
        dof = mesh.edof[ii*mesh.ne+e];

        /* RHS, Division cant be done on element basis because M^-1 ≠ ∑Me^-1 */
        a[dof] = a[dof] + ( f[dof] - C[dof]*vtilde[dof] - prod -
                            KSpring[dof]*dtilde[dof] );

      } /* end dof loop */
    }
  }
 } /* end of element loops */

/* Complete send communications */
MPI_Waitall(8,sendReqs,MPI_STATUS_IGNORE);

#ifdef MY_MPI
mpi.ReqComm(a,false,recReqs,sendReqs,NULL); /* overlapping "forward" communication */
#endif

/* STEP 3: Finish calculating Acceleration and then calculate Displacement and
   Velocity as Correctors */

/* loop INTERNAL dofs by looping INTERNAL nodes */
for(int i=1; i<(mesh.nnx-1); i++){
  for(int j=1; j<(mesh.nny-1); j++){
    for(int k=0; k<(mesh.nnz); k++){

      /* d1 = dof+0, d2 = dof+1, d3 = dof+2 */
      dof = mesh.n(i,j,k)*3;

      /* Divide with global mass vector */
      a[dof+0] /= (M[dof+0] + gamma*dt*C[dof+0]);
      a[dof+1] /= (M[dof+1] + gamma*dt*C[dof+1]);
      a[dof+2] /= (M[dof+2] + gamma*dt*C[dof+2]);

      /* Correct d,v */
      d[dof+0] = dtilde[dof+0];
      d[dof+1] = dtilde[dof+1];
      d[dof+2] = dtilde[dof+2];

      v[dof+0] = vtilde[dof+0] + gamma*dt*a[dof+0];
      v[dof+1] = vtilde[dof+1] + gamma*dt*a[dof+1];
      v[dof+2] = vtilde[dof+2] + gamma*dt*a[dof+2];

      count++;
    }
  }
 }

/* Complete receive communications */
MPI_Waitall(8,recReqs,MPI_STATUS_IGNORE);
mpi.EmptyBuffer(a,false,NULL);

/* loop BOUNDARY dofs by looping BOUNDARY nodes */
for(int i=0; i<(mesh.nnx); i+=(mesh.nnx)){ /* Den er jeg lidt usikker på den her! */
  for(int j=0; j<(mesh.nny); j+=(mesh.nny)){
    for(int k=0; k<(mesh.nnz); k++){

      /* d1 = dof+0, d2 = dof+1, d3 = dof+2 */
      dof = mesh.n(i,j,k)*3;

      /* Divide with global mass vector */
      a[dof+0] /= (M[dof+0] + gamma*dt*C[dof+0]);
      a[dof+1] /= (M[dof+1] + gamma*dt*C[dof+1]);
      a[dof+2] /= (M[dof+2] + gamma*dt*C[dof+2]);

      /* Correct d,v */
      d[dof+0] = dtilde[dof+0];
      d[dof+1] = dtilde[dof+1];
      d[dof+2] = dtilde[dof+2];

      v[dof+0] = vtilde[dof+0] + gamma*dt*a[dof+0];
      v[dof+1] = vtilde[dof+1] + gamma*dt*a[dof+1];
      v[dof+2] = vtilde[dof+2] + gamma*dt*a[dof+2];

      count++;

    }
  }
 }

/* Complete send communications */
MPI_Waitall(8,sendReqs,MPI_STATUS_IGNORE);
