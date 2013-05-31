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

/* For creating and testing for directory */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* #include <vtkStructuredGrid.h> */
/* #include <vtkXMLStructuredGridWriter.h> */
#include <vtkRectilinearGrid.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLPRectilinearGridWriter.h>


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

#include <iostream>
#include <fstream>

void set_structure(FEMclass* mesh,vtkRectilinearGrid *rgrid ){

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

  rgrid->SetDimensions(mesh->nnx,mesh->nny,mesh->nnz);
  rgrid->SetXCoordinates(xCoords);
  rgrid->SetYCoordinates(yCoords);
  rgrid->SetZCoordinates(zCoords);
}

void set_floatVector(FEMclass* mesh,vtkFloatArray *disp,double *vec){

  disp->SetName("disp");
  disp->SetNumberOfComponents(3);
  disp->SetNumberOfTuples(mesh->nn);
  for(int i=0;i<mesh->nn;i++){
    disp->SetTuple3(i,vec[i*3+0], vec[i*3+1], vec[i*3+2]); // set everything to 10
  }

}

#ifdef MY_MPI
void write_file(vtkRectilinearGrid *rgrid, MPIClass *mpi,int step){

    /* Write to file */
  vtkSmartPointer<vtkXMLRectilinearGridWriter>
    writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

  /* create new directory/test for existence */
  /* Consider using the boost library */
  /* http://stackoverflow.com/a/675095/1121523 */
  char dirname [30];
  int n = sprintf(dirname,"plots/sem_%.4i",step);
  struct stat st = {0};
  if (stat(dirname, &st) == -1) {
    mkdir(dirname, 0700);
  }

  char filename [30];
  n = sprintf(filename,"plots/sem_%.4i/sem_%.2i.vtr",step,mpi->rank);
  /* n = sprintf(filename,"plots/sem_%.4i_%i.vtr",step,mpi->rank); */
  writer->SetFileName(filename);

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(rgrid);
#else
  writer->SetInputData(rgrid);
#endif
  writer->Write();


/*   if(mpi->rank==0){ */
/*     vtkSmartPointer<vtkXMLPRectilinearGridWriter> */
/*       pwriter = vtkSmartPointer<vtkXMLPRectilinearGridWriter>::New(); */
/*     char filename2 [30]; */
/*     n = sprintf(filename2,"plots/sem_%.4i.pvtr",step); */

/*     pwriter->SetFileName(filename2); */
/*     pwriter->SetNumberOfPieces(mpi->size); */
/* #if VTK_MAJOR_VERSION <= 5 */
/*     pwriter->SetInput(rgrid); */
/* #else */
/*     pwriter->SetInputData(rgrid); */
/* #endif */
/*     pwriter->Write(); */
/*   } */

}

void write_vtm( MPIClass *mpi,int step){
  /* Write multiblock file */
  /* http://paraview.org/Wiki/ParaView:FAQ#What_file_formats_does_ParaView_support.3F */

  char filename[30];
  int n = sprintf(filename,"plots/sem_%.4i/sem.vtm",step);
  std::ofstream vtkstream(filename);

  vtkstream << "<?xml version=\"1.0\"?>" << "\n";
  vtkstream << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" ";
  vtkstream << "byte_order=\"LittleEndian\" ";
  vtkstream << "compressor=\"vtkZLibDataCompressor\">"  << "\n";
  vtkstream << "	"<<"<vtkMultiBlockDataSet>" << "\n";
  for(int i=0;i<mpi->size;i++){
    n = sprintf(filename,"sem_%.2i.vtr",i);

    vtkstream << "		" <<"<DataSet index=\"" << i << "\"";
    vtkstream << " file=\"" << filename << "\">" << "\n" ;
    vtkstream << "		" <<"</DataSet>" << "\n";
  }

  vtkstream << "	"<<"</vtkMultiBlockDataSet>" << "\n";
  vtkstream << "</VTKFile>" <<"\n";

}

void write_pvd(double dt,int NT){

  char filename[30];
  int n = sprintf(filename,"plots/sem.pvd");
  std::ofstream vtkstream(filename);

  vtkstream << "<?xml version=\"1.0\"?>" << "\n";
  vtkstream << "<VTKFile type=\"Collection\" version=\"0.1\" ";
  vtkstream << "byte_order=\"LittleEndian\" ";
  vtkstream << "compressor=\"vtkZLibDataCompressor\">"  << "\n";
  vtkstream << "	""<Collection>" << "\n";

  for(int i=0;i<NT;i++){
    int n = sprintf(filename,"sem_%.4i/sem.vtm",i);

    vtkstream << "<DataSet timestep=\"" << dt*i << "\" ";
    vtkstream << "file=\"" << filename << "\" />" "\n";
  }
  vtkstream << "	""</Collection>" << "\n";
  vtkstream << "</VTKFile>" <<"\n";
}

void print_vtk(FEMclass* mesh, MPIClass *mpi, int step, double *vec){

  /* Overvej selv at skrive vtu-filerne. */
  /* http://stackoverflow.com/questions/10913666/error-writing-binary-vtk-files */

  /* create new directory/test for existence */
  char dirname [30];
  int n = sprintf(dirname,"plots/");
  struct stat st = {0};
  if (stat(dirname, &st) == -1) {
    mkdir(dirname, 0700);
  }

  vtkSmartPointer<vtkRectilinearGrid>
    rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
  set_structure(mesh,rgrid );

  vtkSmartPointer<vtkFloatArray>
    disp = vtkSmartPointer<vtkFloatArray>::New();
  set_floatVector(mesh,disp,vec);
  /* apply data */
  rgrid->GetPointData()->AddArray(disp);

  write_file(rgrid,mpi,step);
  if(mpi->rank==0){
    write_vtm(mpi,step);
  }


}

#endif



void print_vtk_serial(FEMclass* mesh, int step, double *vec){

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

  //printf("DONE printing\n");
}
