// -*- coding: utf-8 -*-
/* Write the structure to a .vtr file. */
/* VTK (bl. a.Paraview) have two formats. Legacy(old) and XML. Both is a mixture
   of ASCII and binary data. The format used here is rectlinear grid in XML
   format.
   Take a look at vtk-file-formats.pdf for specifications.*/
/* For convince the (huge) VTK library is used. Boyan suggest writing our own
   routines - which is easy enough - but getting to know the VTK library is
   NIIICE!:) */

#include "fedata.h"
#include "vtk.h"

/* #include <vtkStructuredGrid.h> */
/* #include <vtkXMLStructuredGridWriter.h> */
#include <vtkRectilinearGrid.h>
#include <vtkXMLRectilinearGridWriter.h>

/* These are always needed. No matter the mesh type */
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
/* Smart!:) way of initializing variables. Using this when calling, there's no
   need for deleting the pointers manually. Because they are destroyed when this
   routine is finished */
/* http://www.vtk.org/Wiki/VTK/Tutorials/SmartPointers */
#include <vtkSmartPointer.h>
#include <vtkPointData.h> /* Needed for adding data to the mesh */



void print_vtk(FEMclass* mesh, int step, double *vec){

  /* USE SMARTPOINTERS */

  // Create a rectilinear grid by defining three arrays specifying the
  // coordinates in the x-y-z directions.
  vtkSmartPointer<vtkDoubleArray>
    xCoords = vtkSmartPointer<vtkDoubleArray>::New();
  for (int i=0; i<mesh->nnx; i++) xCoords->InsertNextValue(mesh->x[i]);
  vtkSmartPointer<vtkDoubleArray>
    yCoords = vtkSmartPointer<vtkDoubleArray>::New();
  for (int i=0; i<mesh->nny; i++) yCoords->InsertNextValue(mesh->y[i]);
  vtkSmartPointer<vtkDoubleArray>
    zCoords = vtkSmartPointer<vtkDoubleArray>::New();
  for (int i=0; i<mesh->nnz; i++) zCoords->InsertNextValue(mesh->z[i]);


  // The coordinates are assigned to the rectilinear grid. Make sure that
  // the number of values in each of the XCoordinates, YCoordinates,
  // and ZCoordinates is equal to what is defined in SetDimensions().
  vtkSmartPointer<vtkRectilinearGrid>
    rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
  rgrid->SetDimensions(mesh->nnx,mesh->nny,mesh->nnz);
  rgrid->SetXCoordinates(xCoords);
  rgrid->SetYCoordinates(yCoords);
  rgrid->SetZCoordinates(zCoords);


  /* Write data to the grid */
  /* vtkFloatArray* temperature = vtkFloatArray::New(); */
  /* vtkSmartPointer<vtkFloatArray> */
  /*   temperature = vtkSmartPointer<vtkFloatArray>::New(); */
  /* temperature->SetName("Temperature"); */
  /* temperature->SetNumberOfComponents(1); */
  /* temperature->SetNumberOfValues(mesh.nn); */
  /* for(int i=0;i<mesh.nn;i++){ */
  /*   temperature->SetValue(i,i); */
  /* } */
  /* rgrid->GetPointData()->Addrray(temperature); */

  /* vtkSmartPointer<vtkFloatArray> */
  /*   velocity = vtkSmartPointer<vtkFloatArray>::New(); */
  /* velocity->SetName("Velocity"); */
  /* velocity->SetNumberOfComponents(3); */
  /* velocity->SetNumberOfTuples(mesh.nn); */
  /* for(int i=0;i<mesh.nn;i++){ */
  /*   velocity->SetTuple3(i,1*i, 10, 30); // set everything to 10 */
  /* } */
  /* rgrid->GetPointData()->AddArray(velocity); */

  vtkSmartPointer<vtkFloatArray>
    disp = vtkSmartPointer<vtkFloatArray>::New();
  disp->SetName("disp");
  disp->SetNumberOfComponents(3);
  disp->SetNumberOfTuples(mesh->nn);
  for(int i=0;i<mesh->nn;i++){
    disp->SetTuple3(i,vec[i*3+0], vec[i*3+1], vec[i*3+2]); // set everything to 10
  }
  rgrid->GetPointData()->AddArray(disp);


  /* Write to file */
  vtkSmartPointer<vtkXMLRectilinearGridWriter>
    writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

  char filename [30];
  int n = sprintf(filename,"plots/sem_%.4i.vtr",step);
  writer->SetFileName(filename);

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(rgrid);
#else
  writer->SetInputData(rgrid);
#endif
  writer->Write();

  /* clean up */
  /* xCoords->Delete(); */
  /* yCoords->Delete(); */
  /* zCoords->Delete(); */
  /* rgrid->Delete(); */

  printf("DONE printing\n");
}
