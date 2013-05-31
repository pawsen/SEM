/* STEP 1: Define Predictors */

for(int i=0; i<(mesh.nn)*3; i++){

  count++;

  dtilde[i] = d[i] + v[i]*dt + 0.5*dt*dt*a[i];
  vtilde[i] = v[i] + (1-gamma)*dt*a[i];

  a[i] = 0.0; // while looping, also reset acceleration for step 2

 }

/* "backward" communication */
#ifdef MY_MPI
mpi.communicate(dtilde,true);
mpi.communicate(vtilde,true);
#endif


/* STEP 2: Solve for Acceleration */

/* loop elements, e */
int dof; // single entry of edof-vector
double prod = 0; // single entry of ke*dtilde product

for(int e=0; e<mesh.ne; e++){

  /* loop dofs, i in edofs-vec for this element */
  for(int i=0; i<mesh.nen*3; i++){

    /* prod = dof-th entry of vector from ke*dtilde product */
    prod = 0;
    for(int j=0; j<mesh.nen*3; j++){
      count++;
      prod = prod + Ke[i*mesh.nen*3+j]*dtilde[ mesh.edof[j*mesh.ne+e] ];
    }
    // global index of the considered dof of the element
    dof = mesh.edof[i*mesh.ne+e];

    /* RHS, Division cant be done on element basis because M^-1 ≠ ∑Me^-1 */
    a[dof] = a[dof] + ( f[dof] - C[dof]*vtilde[dof] - prod -
                        KSpring[dof]*dtilde[dof] );

  } /* end dof loop */
 } /* end of element loops */

#ifdef MY_MPI
mpi.communicate(a,false); /* "forward" communication */
#endif

/* Divide with global mass vector */
for(int i=0;i<mesh.nn*3;i++){
  count++;
  a[i] /= (M[i] + gamma*dt*C[i]);
 }

/* STEP 3: Calculate Displacement and Velocity as Correctors */

for(int i=0; i<mesh.nn*3; i++){
  count++;
  d[i] = dtilde[i];
  v[i] = vtilde[i] + gamma*dt*a[i];
 }
