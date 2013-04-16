// -*- coding: utf-8 -*-

#include "mpidata.h"
#include <math.h>
#include <iostream>
using namespace std;

MPIClass::MPIClass (int argc, char **argv,DataStruct *data_in){

  data = data_in;
  initMPI(argc,argv);
  partitioning();
  getNeigh();

}

void MPIClass::initMPI(int argc, char **argv){
  /**************************/
  /* STEP 0: Initialize MPI */
  /**************************/

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); //my rank
  MPI_Comm_size(MPI_COMM_WORLD, &size); //proc size

}

void MPIClass::partitioning(){

  /*****************/
  /* Partitioning  */
  /*****************/

  int gnelx=data->inInt[0],gnely=data->inInt[1],gnelz=data->inInt[2];
  double glx=data->inDouble[0],gly=data->inDouble[1],glz=data->inDouble[2];


  /* 2D partitioning */
  const int ndims = 2;
  /* Is the grid periodic. 0: no, 1: yes */
  int periods[ndims] = {0,0};
  int nproc = size;
  /* Dims_create generates the decomposition. Eg. a call with nproc = 6 returns
     dims={3,2}. This correspond to a 2x3 grid. MPI uses c-style row-major
     numbering, but the first entry in dims correspond to the number of columns.
     STUPID, but that's life. With dims={3,2} Cart_create returns following
     coordinates and rank:

     (0,0)=0 (0,1)=1 (0,2)=2
     (1,0)=3 (1,1)=4 (1,2)=5
     Info: http://stackoverflow.com/a/13783336/1121523
  */
  /* Let MPI decide the distribution(dims = {0}). If we want a fixed number of
     dims in one direction, it can be specified here*/
  dims[0] = 0; dims[1] = 0;
  MPI_Dims_create(nproc,ndims,dims);
  int nProcx=dims[1], nProcy=dims[0];

  /* Create cartesian communicator */
  MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,0,&comm_cart);
  if (rank == 0)
    cout << "rank: " << rank << ", dims: " << dims[0] << " " << dims[1] << endl;

  /* Local dimensions of mesh */
  int dny = gnely % nProcy;   /* remainder of mesh partitioning along x-axis */
  int dnx = gnelx % nProcx;   /* remainder of mesh partitioning along y-axis */
  int nelx = ceil(gnelx/nProcx);
  int nely = ceil(gnely/nProcy);
  int nelz = gnelz;

  /* Determine cartesian coords of this process */
  MPI_Cart_coords(comm_cart,rank,ndims,myCoords);
  cout << "rank: " << rank << ", coor: " << myCoords[0] << " " << myCoords[1] << endl;


  /* Distribute remainder elements on to first dnx,dny domain compositions */
  if( myCoords[0] < dny )
    nely++;

  if( myCoords[1] < dnx )
    nelx++;

  /* Physical dimensions of local mesh. Total length */
  int lx = glx/nProcx;
  int ly = gly/nProcy;
  data->outDouble[0] = lx;/* lx */
  data->outDouble[1] = ly; /* ly */
  data->outDouble[2] = glz;        /* lz */
  data->outInt[0] = nelx;
  data->outInt[1] = nely;
  data->outInt[2] = nelz;

  /* Mesh offset in x/y direction */
  offsetX = myCoords[1] * lx;
  offsetY = myCoords[0] * ly;
  /* cout << "myCoorX: "  << myCoords[1] << ", offsetX: " << offsetX */
  /*      << ", myCoorY: "  << myCoords[0] << ", offsetY: " << offsetY << endl; */
  /* if (myCoords[0] < dny) */
  /*   offsetY = myCoords[0]* (nely+1); */
  /* else */
  /*   offsetY = myCoords[0]*nely+dny; */

  /* if (myCoords[1] < dnx) */
  /*   offsetX = myCoords[1]* (nelx+1); */
  /* else */
  /*   offsetX = myCoords[1]*nelx+dnx; */

  if(rank == 0 ){
    cout << "rank: " << rank << ", lx,ly,lz: ";
    for(int i=0;i<3;i++){
      cout << data->outDouble[i] << " ";
    }
    cout << endl << "rank: " << rank << ", nelx,nely,nelz: ";
    for(int i=0;i<3;i++){
      cout << data->outInt[i] << " ";
    }
    cout << endl;
  }
}


void MPIClass::getNeigh(){

  /* Determine neighbour procs.
     Each proc, must send to "dest" to and receive from "src".
     Out-of-range procs will result in MPI_PROC_NULL when shifting
     outside domain and not using periodicity.
  */
  /*
    Sending to MPI_PROC_NULL(integer = -2) result in:
    A communication with process MPI_PROC_NULL has no effect. A send to
    MPI_PROC_NULL succeeds and returns as soon as possible. A receive from
    MPI_PROC_NULL succeeds and returns as soon as possible with no
    modifications to the receive buffer
  */
  /*
    So for nodes on the boundary, we need to determine the Out-of-range
    coordinates manually, since MPI_Cart_rank return:

    Out-of-range coordinates are erroneous for non-periodic dimensions. Versions
    of MPICH before 1.2.2 returned MPI_PROC_NULL for the rank in this case.
    Newer versions abort execution.

  */

  int tmpCoord[2];
  /* {N,NE,E,SE,S, SW,W,NW} */
  int tmpCoord1[8] = {-1,-1,0,1, 1, 1, 0,-1};
  int tmpCoord2[8] = { 0, 1,1,1, 0,-1,-1,-1};

  for(int i=0;i<8;i++){
    tmpCoord[0] = myCoords[0] + tmpCoord1[i];
    tmpCoord[1] = myCoords[1] + tmpCoord2[i];

    if (rank == 0){
      // cout << tmpCoord[i] << " " << tmpCoord[1] << " | ";
    }

    if(tmpCoord[0]<0){ neigh[i]=MPI_PROC_NULL; }
    else if(tmpCoord[0]>(dims[0]-1)){ neigh[i]=MPI_PROC_NULL; }
    else if(tmpCoord[1]<0){ neigh[i]=MPI_PROC_NULL; }
    else if(tmpCoord[1]>(dims[1]-1)){ neigh[i]=MPI_PROC_NULL; }
    else {
      MPI_Cart_rank(comm_cart,tmpCoord,&neigh[i]);
    }
  }
  /* Works! */
  /* MPI_Barrier; */
  /* fflush(stdout); */
  if (rank < 8){
    cout << "rank: "<< rank << ", neigh: ";
    for(int i=0;i<8;i++){
      cout << neigh[i] << ", ";
    }
    cout << endl;
  }

}

void MPIClass::datatype(){


  /**************************************************/
  /* Get displacement of each arrays starting point */
  /**************************************************/
  /* http://stackoverflow.com/a/2655363/1121523 */

  int N=10;            // number of arrays (first dimension)
  int sizes[N];     // number of elements in each array (second dimensions)
  int* arrays[N];   // pointers to the start of each array
  /* Get starting point */
  MPI_Aint base;
  MPI_Address(arrays[0], &base);

  /* Get displacement */
  MPI_Aint* displacements = new MPI_Aint[N];
  for (int i=0; i<N; ++i){
    /* Brug MPI_Get_address i stedet.
       http://mpi.deino.net/mpi_functions/MPI_Get_address.html */
    MPI_Address(arrays[i], &displacements[i]);
    displacements[i] -= base;
  }

  /* construct data type */
  MPI_Datatype newType;
  MPI_Type_hindexed(N, sizes, displacements, MPI_INTEGER, &newType);
  MPI_Type_commit(&newType);


}

void MPIClass::datatype_helper(){
  int node, stride, nblock,ntot,iNode,jNode,kNode;
  int nnx = mesh->nnx, nny = mesh->nny, nnz = mesh->nnz;

  /* These two are the in/out of plane. And are not communicated */
  /* Z-top */
  /* DOF =  nnx*nny*(nnz-1)*3..nnx*nny*nnz*3-1 */
  /* Z-bottom */
  /* DOF =  0..(nnx*nny)*3 -1*/

  int orientation = N;
  switch(orientation){
  case N:

    /* North - xz plane for y=nny-1 */
    jNode = nny -1;    /* REMEMBER! zero based */
    /* All x-dofs */
    for(int i=0;i<nnx-1;i++){
      for(int k=0;k<nnz-1;k++){
        node = mesh->n(i,jNode,k)*3;
      }
    }
    /* first x-dof */
    node = mesh->n(0,jNode,0)*3;
    /* stride,  eg. length of continuous block */
    stride = nnx*3;
    /* number of blocks */
    nblock = nnz;
    /* total number of DOF's to be send */
    ntot = stride*nblock;
    break;
  case E:

    break;
  }
}

void MPIClass::initBuffers(FEMclass* mesh_in){

  mesh = mesh_in;

  /* init send/receive buffers */
  int nnx = mesh->nnx, nny = mesh->nny, nnz = mesh->nnz;

  //buf = new double*[8];

  lengthBuf[N] = nnx*nnz*3;
  buf[N] = new double[lengthBuf[N]];

  lengthBuf[NE] = nnz*3;
  buf[NE] = new double[lengthBuf[NE]];

  lengthBuf[E] = nny*nnz*3;
  buf[E] = new double[lengthBuf[E]];

  lengthBuf[SE] = nnz*3;
  buf[SE] = new double[lengthBuf[SE]];

  lengthBuf[S] = nnx*nnz*3;
  buf[S] = new double[lengthBuf[S]];

  lengthBuf[SW] = nnz*3;
  buf[SW] = new double[lengthBuf[SW]];

  lengthBuf[W] = nny*nnz*3;
  buf[W] = new double[lengthBuf[W]];

  lengthBuf[NW] = nnz*3;
  buf[NW] = new double[lengthBuf[NW]];

}

void MPIClass::fillBuffer(double *data, int face){

  //cout << "rank:" << rank << ", fill buffer for face:" << face << endl;

  int nnx = mesh->nnx, nny = mesh->nny, nnz = mesh->nnz;
  int iNode,jNode,kNode;
  int ii;

  // cout << "rank:" << rank << ", fill face: " << face << endl;

  switch(face){

  case N:
    /* North: (i,j,k) = (i,0,k) */
    jNode = 0;
    ii=0;
    for(int i=0; i<mesh->nnx; i++){
      for(int k=0; k<mesh->nnz; k++){
        buf[N][ii+0] = data[mesh->n(i,jNode,k)*3 + 0 ];
        buf[N][ii+1] = data[mesh->n(i,jNode,k)*3 + 1 ];
        buf[N][ii+2] = data[mesh->n(i,jNode,k)*3 + 2 ];
        ii += 3;
      }
    }
    break;

  case NE:
    /* North East: (i,j,k) = (nnx-1,0,k) */
    iNode = mesh->nnx-1;
    jNode = 0;
    ii=0;
    for(int k=0; k<mesh->nnz; k++){
      buf[NE][ii+0] = data[mesh->n(iNode,jNode,k)*3 + 0 ];
      buf[NE][ii+1] = data[mesh->n(iNode,jNode,k)*3 + 1 ];
      buf[NE][ii+2] = data[mesh->n(iNode,jNode,k)*3 + 2 ];
      ii += 3;
    }
    break;

  case E:
    /* East (i,j,k) = (nnx-1,j,k) */
    iNode = mesh->nnx-1;
    ii=0;
    for(int j=0; j<mesh->nny; j++){
      for(int k=0; k<mesh->nnz; k++){
        buf[E][ii+0] = data[mesh->n(iNode,j,k)*3 + 0 ];
        buf[E][ii+1] = data[mesh->n(iNode,j,k)*3 + 1 ];
        buf[E][ii+2] = data[mesh->n(iNode,j,k)*3 + 2 ];
        ii += 3;
      }
    }
    break;

  case SE:
    /* South East: (i,j,k) = (nnx-1,nny-1,k) */
    iNode = mesh->nnx-1;
    jNode = mesh->nny-1;
    ii=0;
    for(int k=0; k<mesh->nnz; k++){
      buf[SE][ii+0] = data[mesh->n(iNode,jNode,k)*3 + 0 ];
      buf[SE][ii+1] = data[mesh->n(iNode,jNode,k)*3 + 1 ];
      buf[SE][ii+2] = data[mesh->n(iNode,jNode,k)*3 + 2 ];
      ii += 3;
    }
    break;

  case S:
    /* South: (i,j,k) = (i,nny-1,k) */
    jNode = mesh->nny-1;
    ii=0;
    for(int i=0; i<mesh->nnx; i++){
      for(int k=0; k<mesh->nnz; k++){
        buf[S][ii+0] = data[mesh->n(i,jNode,k)*3 + 0 ];
        buf[S][ii+1] = data[mesh->n(i,jNode,k)*3 + 1 ];
        buf[S][ii+2] = data[mesh->n(i,jNode,k)*3 + 2 ];
        ii += 3;
      }
    }
    break;

  case SW:
    /* South West: (i,j,k) = (0,nny-1,k) */
    iNode = 0;
    jNode = mesh->nny-1;
    ii=0;
    for(int k=0; k<mesh->nnz; k++){
      buf[SW][ii+0] = data[mesh->n(iNode,jNode,k)*3 + 0 ];
      buf[SW][ii+1] = data[mesh->n(iNode,jNode,k)*3 + 1 ];
      buf[SW][ii+2] = data[mesh->n(iNode,jNode,k)*3 + 2 ];
      ii += 3;
    }
    break;

  case W:
    /* West: (i,j,k) = (0,j,k) */
    iNode = 0;
    ii=0;
    for(int j=0; j<mesh->nny; j++){
      for(int k=0; k<mesh->nnz; k++){
        buf[W][ii+0] = data[mesh->n(iNode,j,k)*3 + 0 ];
        buf[W][ii+1] = data[mesh->n(iNode,j,k)*3 + 1 ];
        buf[W][ii+2] = data[mesh->n(iNode,j,k)*3 + 2 ];
        ii += 3;
      }
    }
    break;

  case NW:
    /* North West: (i,j,k) = (0,0,k) */
    iNode = 0;
    jNode = 0;
    ii=0;
    for(int k=0; k<mesh->nnz; k++){
      buf[NW][ii+0] = data[mesh->n(iNode,jNode,k)*3 + 0 ];
      buf[NW][ii+1] = data[mesh->n(iNode,jNode,k)*3 + 1 ];
      buf[NW][ii+2] = data[mesh->n(iNode,jNode,k)*3 + 2 ];
      ii += 3;
    }
    break;

  } /* end switch */

}


void MPIClass::addBuffer(double *data, int face){

  //cout << "rank:" << rank << ", add buffer to face:" << face << endl;

  int nnx = mesh->nnx, nny = mesh->nny, nnz = mesh->nnz;
  int iNode,jNode,kNode;
  int ii;

  switch(face){

  case N:
    /* North: (i,j,k) = (i,0,k) */
    jNode = 0;
    ii=0;
    for(int i=0; i<mesh->nnx; i++){
      for(int k=0; k<mesh->nnz; k++){
        data[mesh->n(i,jNode,k)*3 + 0 ] += buf[N][ii+0];
        data[mesh->n(i,jNode,k)*3 + 1 ] += buf[N][ii+1];
        data[mesh->n(i,jNode,k)*3 + 2 ] += buf[N][ii+2];
        ii += 3;
      }
    }
    break;

  case NE:
    /* North East: (i,j,k) = (nnx-1,0,k) */
    iNode = mesh->nnx-1;
    jNode = 0;
    ii=0;
    for(int k=0; k<mesh->nnz; k++){
      data[mesh->n(iNode,jNode,k)*3 + 0 ] += buf[NE][ii+0];
      data[mesh->n(iNode,jNode,k)*3 + 1 ] += buf[NE][ii+1];
      data[mesh->n(iNode,jNode,k)*3 + 2 ] += buf[NE][ii+2];
      ii += 3;
    }
    break;

  case E:
    /* East (i,j,k) = (nnx-1,j,k) */
    iNode = mesh->nnx-1;
    ii=0;
    for(int j=0; j<mesh->nny; j++){
      for(int k=0; k<mesh->nnz; k++){
        data[mesh->n(iNode,j,k)*3 + 0 ] += buf[E][ii+0];
        data[mesh->n(iNode,j,k)*3 + 1 ] += buf[E][ii+1];
        data[mesh->n(iNode,j,k)*3 + 2 ] += buf[E][ii+2];
        ii += 3;
      }
    }
    break;

  case SE:
    /* South East: (i,j,k) = (nnx-1,nny-1,k) */
    iNode = mesh->nnx-1;
    jNode = mesh->nny-1;
    ii=0;
    for(int k=0; k<mesh->nnz; k++){
      data[mesh->n(iNode,jNode,k)*3 + 0 ] += buf[SE][ii+0];
      data[mesh->n(iNode,jNode,k)*3 + 1 ] += buf[SE][ii+1];
      data[mesh->n(iNode,jNode,k)*3 + 2 ] += buf[SE][ii+2];
      ii += 3;
    }
    break;

  case S:
    /* South: (i,j,k) = (i,nny-1,k) */
    jNode = mesh->nny-1;
    ii=0;
    for(int i=0; i<mesh->nnx; i++){
      for(int k=0; k<mesh->nnz; k++){
        data[mesh->n(i,jNode,k)*3 + 0 ] += buf[S][ii+0];
        data[mesh->n(i,jNode,k)*3 + 1 ] += buf[S][ii+1];
        data[mesh->n(i,jNode,k)*3 + 2 ] += buf[S][ii+2];
        ii += 3;
      }
    }
    break;

  case SW:
    /* South West: (i,j,k) = (0,nny-1,k) */
    iNode = 0;
    jNode = mesh->nny-1;
    ii=0;
    for(int k=0; k<mesh->nnz; k++){
      data[mesh->n(iNode,jNode,k)*3 + 0 ] += buf[SW][ii+0];
      data[mesh->n(iNode,jNode,k)*3 + 1 ] += buf[SW][ii+1];
      data[mesh->n(iNode,jNode,k)*3 + 2 ] += buf[SW][ii+2];
      ii += 3;
    }
    break;

  case W:
    /* West: (i,j,k) = (0,j,k) */
    iNode = 0;
    ii=0;
    for(int j=0; j<mesh->nny; j++){
      for(int k=0; k<mesh->nnz; k++){
        data[mesh->n(iNode,j,k)*3 + 0 ] += buf[W][ii+0];
        data[mesh->n(iNode,j,k)*3 + 1 ] += buf[W][ii+1];
        data[mesh->n(iNode,j,k)*3 + 2 ] += buf[W][ii+2];
        ii += 3;
      }
    }
    break;

  case NW:
    /* North West: (i,j,k) = (0,0,k) */
    iNode = 0;
    jNode = 0;
    ii=0;
    for(int k=0; k<mesh->nnz; k++){
      data[mesh->n(iNode,jNode,k)*3 + 0 ] += buf[NW][ii+0];
      data[mesh->n(iNode,jNode,k)*3 + 1 ] += buf[NW][ii+1];
      data[mesh->n(iNode,jNode,k)*3 + 2 ] += buf[NW][ii+2];
      ii += 3;
    }
    break;

  } /* end switch */



}


void MPIClass::overwriteBuffer(double *data, int face){

  /****************************************************************************/
  /* A process should only receive from processes with higher ranks, which owns
     the boundary nodes, when having nodal values overwritten at boundarie s. */
  /* Therefore it is only possible to to receive from E, SE, S neighbours.    */
  /* If: no neighbours are apparent to SE and S -> E neighbour owns all the
     boundary nodes nodes on E plane. */
  /* Else if: no neighbour is apparent to SE -> S neighbour owns all boundary
     nodes on S plane. */
  /* Else: SE is apparent and it owns the nodes on the SE corner edge, why these
     should not be overwritten by S. */
  /***************************************************************************/

  //  cout << "rank:" << rank << ", overwrite with buffer on face:" << face << endl;

  int nnx = mesh->nnx, nny = mesh->nny, nnz = mesh->nnz;
  int iNode,jNode,kNode;
  int iMax,jMax;
  int ii;

  switch(face){

  case E:

    /* check for SE or S neighbour */
    if( neigh[SE]!=MPI_PROC_NULL && neigh[S]!=MPI_PROC_NULL )
      jMax = mesh->nny-1;   /* exclude corner */
    else
      jMax = mesh->nny;

    /* East (i,j,k) = (nnx-1,j,k) */
    iNode = mesh->nnx-1;
    ii=0;
    for(int j=0; j<jMax; j++){
      for(int k=0; k<mesh->nnz; k++){
        data[mesh->n(iNode,j,k)*3 + 0 ] = buf[E][ii+0];
        data[mesh->n(iNode,j,k)*3 + 1 ] = buf[E][ii+1];
        data[mesh->n(iNode,j,k)*3 + 2 ] = buf[E][ii+2];
        ii += 3;
      }
    }
    break;

  case SE:
    /* South East: (i,j,k) = (nnx-1,nny-1,k) */
    iNode = mesh->nnx-1;
    jNode = mesh->nny-1;
    ii=0;
    for(int k=0; k<mesh->nnz; k++){
      data[mesh->n(iNode,jNode,k)*3 + 0 ] = buf[SE][ii+0];
      data[mesh->n(iNode,jNode,k)*3 + 1 ] = buf[SE][ii+1];
      data[mesh->n(iNode,jNode,k)*3 + 2 ] = buf[SE][ii+2];
      ii += 3;
    }
    break;

  case S:

    /* /\* check for SE neighbour *\/ */
    if( neigh[SE]!=MPI_PROC_NULL )
      iMax = mesh->nnx-1;   /* exclude corner */
    else
      iMax = mesh->nnx;

    /* South: (i,j,k) = (i,nny-1,k) */
    jNode = mesh->nny-1;
    ii=0;
    for(int i=0; i<iMax; i++){
      for(int k=0; k<mesh->nnz; k++){
        data[mesh->n(i,jNode,k)*3 + 0 ] = buf[S][ii+0];
        data[mesh->n(i,jNode,k)*3 + 1 ] = buf[S][ii+1];
        data[mesh->n(i,jNode,k)*3 + 2 ] = buf[S][ii+2];
        ii += 3;
      }
    }
    break;

  } /* end switch */

}



void MPIClass::communicate(double *data, bool backward){

  MPI_Request recReqs[8], sendReqs[8];
  int sendTag = 0;
  int numRecReqs = 0, numSendReqs = 0;

  for(int face=0; face<8; face++){

    /* Default should be null request */
    recReqs[face] = MPI_REQUEST_NULL;
    sendReqs[face] = MPI_REQUEST_NULL;

    if (neigh[face] == MPI_PROC_NULL); /* do nothing */
    else if( (neigh[face] < rank && !backward) || (neigh[face] > rank && backward) ){ // request a receive

      numRecReqs++;

      MPI_Irecv(&buf[face][0],lengthBuf[face],MPI_DOUBLE,neigh[face],MPI_ANY_TAG,comm_cart,&recReqs[face]);

    } /* end else if */

    else{ // request a send

      numSendReqs++;

      fillBuffer(data,face);    /* fill only buffers for relevant faces */

      MPI_Isend(&buf[face][0],lengthBuf[face],MPI_DOUBLE,neigh[face],sendTag,comm_cart,&sendReqs[face]);

    } /*  end else */

  }   /* end for */

  // Receive data from other procs
  int i, numRec = 0;
  while(numRec<numRecReqs){

    numRec++; // increment number of receipts completede

    MPI_Waitany(8,recReqs,&i,MPI_STATUS_IGNORE); /* i is index of the completed receipt */

    if(backward)
      overwriteBuffer(data,i); /* overwrite vector by received contributions while waiting */
    else
      addBuffer(data,i); /* add received contributions to vector while waiting */

  }

}
