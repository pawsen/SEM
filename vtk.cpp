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

#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkRectilinearGrid.h>
#include <vtkXMLRectilinearGridWriter.h>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include "vtkPointData.h"

//#include <vtkXMLUnstructuredGridWriter.h>


void print_vtk(FEMclass mesh, int step){

  /* USE SMARTPOINTERS */
  
  /* http://www.vtk.org/Wiki/VTK/Tutorials/SmartPointers */

  // Create a rectilinear grid by defining three arrays specifying the
  // coordinates in the x-y-z directions.
  vtkDoubleArray *xCoords = vtkDoubleArray::New();
  for (int i=0; i<mesh.nnx; i++) xCoords->InsertNextValue(i);
 
  vtkDoubleArray *yCoords = vtkDoubleArray::New();
  for (int i=0; i<mesh.nny; i++) yCoords->InsertNextValue(i);
  
  vtkDoubleArray *zCoords = vtkDoubleArray::New();
  for (int i=0; i<mesh.nnz; i++) zCoords->InsertNextValue(i);


  // The coordinates are assigned to the rectilinear grid. Make sure that
  // the number of values in each of the XCoordinates, YCoordinates,
  // and ZCoordinates is equal to what is defined in SetDimensions().
  vtkRectilinearGrid *rgrid = vtkRectilinearGrid::New();
  rgrid->SetDimensions(mesh.nnx,mesh.nny,mesh.nnz);
  rgrid->SetXCoordinates(xCoords);
  rgrid->SetYCoordinates(yCoords);
  rgrid->SetZCoordinates(zCoords);


  /* Write data to the grid */
  vtkFloatArray* temperature = vtkFloatArray::New();
  temperature->SetName("Temperature");
  temperature->SetNumberOfComponents(1);
  temperature->SetNumberOfValues(mesh.nn);
  for(int i=0;i<mesh.nn;i++){
    temperature->SetValue(i,i);
  }
  rgrid->GetPointData()->AddArray(temperature);

  vtkFloatArray* velocity = vtkFloatArray::New();
  velocity->SetName("Velocity");
  velocity->SetNumberOfComponents(3);
  velocity->SetNumberOfTuples(mesh.nn);
  for(int i=0;i<mesh.nn;i++){
    velocity->SetTuple3(i,1*i, 10, 30); // set everything to 10
  }
  rgrid->GetPointData()->AddArray(velocity);


  // Create a double array.
  /* vtkDoubleArray* vorticity = vtkDoubleArray::New(); */
  /* vorticity->SetName("Vorticity"); */
  /* vorticity->SetNumberOfComponents(3); */
  /* vorticity->SetNumberOfTuples(mesh.nn*3); */
  /* vorticity->InsertNextValue(3.4); */

  // Create the dataset. In this case, we create a vtkPolyData
  /* vtkPolyData* polydata = vtkPolyData::New(); */
  /* polydata->GetPointData()->SetScalars(temperature); */
  // Assign scalars
  /* rgrid->GetPointData()->SetScalars(temperature); */
  // Add the vorticity array. In this example, this field
  // is not used.
  // rgrid->GetPointData()->AddArray(vorticity);

  /* Write to file */
  vtkSmartPointer<vtkXMLRectilinearGridWriter>
    writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

  char filename [30];
  int n = sprintf(filename,"sem_%.4i.vtr",step);
  writer->SetFileName(filename);

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(rgrid);
#else
  writer->SetInputData(rgrid);
#endif
  writer->Write();

  /* clean up */
  xCoords->Delete();
  yCoords->Delete();
  zCoords->Delete();
  rgrid->Delete();

  printf("DONE printing\n");
}
