#ifndef VTK_H
#define VTK_H

#ifdef MY_MPI
#include "mpidata.h"
#endif


void print_vtk_serial(FEMclass *mesh, int step, double *vect);

#ifdef MY_MPI
void write_pvd(double dt,int NT);
void print_vtk(FEMclass* mesh, MPIClass *mpi, int step, double *vec);
#endif



#endif
