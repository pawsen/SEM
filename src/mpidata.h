// -*- coding: utf-8 -*-

#ifndef MPIDATA_H
#define MPIDATA_H

#include "fedata.h"
#include <mpi.h>


/* Why use typdedef: http://stackoverflow.com/a/1675446/1121523 */
typedef struct{
  int inInt[3], outInt[3];
  double inDouble[3], outDouble[3];
}DataStruct;

class MPIClass{
private:

  FEMclass *mesh;
  DataStruct *data;

  enum {N=0,NE,E,SE,S, SW,W,NW};

  void initMPI(int argc, char **argv);
  void partitioning();
  /* Get neighbors to current proc */
  void getNeigh();
  void fillBuffer();
  void datatype_construct(int firstDOF, int stride,int distBlocks,
                          int nblocks, double *data);

  void datatype_helper(MPI_Datatype newType,int face, /*in */
                       int firstDOF, int stride, int distBlocks,
                       int nblocks );

  /* Array holding the proc id of neighbors to current proc.
     If no neighbor the value is MPI_PROC_NULL==-2*/
  int neigh[8];

  /* Sending/receiving buffers using a array-of-array decleration allowing more compact notation */
  /* Array of the 8 pointers which points toward Nbuf, NEbuf... see initBuffer */
  double *buf[8];
  /* Lenghts of the 8 buffers */
  int lengthBuf[8];

public:

  /* MPI data types and related */
  /* Cartesian communicator */
  MPI_Comm comm_cart;

  /* Constructor */
  MPIClass (int argc, char **argv,DataStruct *data_in);

  // void getComm(FEMclass *mesh)

  /* Init the send/recv buffer by letting the  MPIclass see the mesh */
  void initBuffers(FEMclass*);

  /* Fill buffer with data to send */
  void fillBuffer(double*,int,int);

  /* Add from buffer with received data */
  void addBuffer(double*,int,int);

  /* Overwrite data by data in buffer */
  void overwriteBuffer(double *,int,int);

  /* Communicate nodal values at boundaries to neighbours */
  void communicate(double *,bool);

  /* Overlapping communications */
  /* fill buffers and requests sends+recs */
  void ReqComm(double *,bool,MPI_Request[8],MPI_Request[8],double*);
  /* to be used after ReqComm */
  void EmptyBuffer(double *,bool,double*);

  /* Coordinates of local proc  */
  int myCoords[2];

  int dims[2];
  int rank, size;
  /* Mesh offset in x/y direction */
  int offsetX, offsetY;
};

//extern DataStruct data;

#endif
