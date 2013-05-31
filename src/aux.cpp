// -*- coding: utf-8 -*-

#include "aux.h"
#include "fedata.h"
#include "libelement.h"

/* Use 'extern "C"' to tell the C++ compiler that the included library is
   compiled with C and not C++
   http://stackoverflow.com/a/67930/1121523 */
extern "C" {
#include <cblas.h>
  //#include <atlas/cblas.h>
}

#include <iostream>
using namespace std;


void matmul_test(){

  double *A = new double[4*3];
  double *B = new double[4*5];
  double *C = new double[3*5];

  int *shape_A = new int[2];
  int *shape_B = new int[2];
  int *shape_C = new int[2];

  /*************/
  /* TRANSPOSE */
  /*************/
  A[0] = 14.0; A[1] = 9.0;  A[2] = 3.0;     // 4x3
  A[3] = 2.0;  A[4] = 11.0; A[5] = 15.0;
  A[6] = 0.0;  A[7] = 12.0; A[8] = 17.0;
  A[9] = 5.0;  A[10]= 2.0;  A[11]= 3.0;

  B[0] = 3; B[1]=  7; B[2] =12; B[3] =22; B[4] =1;
  B[5] =45; B[6] = 2; B[7] =19; B[8] =10; B[9] =5;
  B[10]= 8; B[11]=23; B[12]=19; B[13]= 4; B[14]=7;
  B[15]= 9; B[16]=16; B[17]= 3 ;B[18]=10; B[19]=2;

  /* Remember: The dimensions(shape_) is the dimension of the matrix used in the
     multiplication. So shape_A for A and A' are different*/
  shape_A[0] = 3; shape_A[1] = 4;
  shape_B[0] = 4; shape_B[1] = 5;
  shape_C[0] = shape_A[0]; shape_C[1] = shape_B[1];
  matmul(1.0,A,B,0.0,C,shape_A,shape_B,true);

  /* C = */
  /* 177   182   221   378    34 */
  /* 636   393   551   376   152 */
  /* 847   490   653   314   203 */

  cout << "aT*b=C-matric from blas: \n";
  for(int i=0;i<shape_C[0];i++){
    for(int j=0;j<shape_C[1];j++){
      cout << C[i*shape_C[1]+j] << "\t";
    }
    cout << endl;
  }

  delete[] B; delete[] C;

  B = new double[6];
  C = new double[8];

  /**********/
  /* NORMAL */
  /**********/
  B[0] = 12.0;  B[1] = 25.0;                // 3x2
  B[2] = 9.0;   B[3] = 10.0;
  B[4] = 8.0;   B[5] = 5.0;

  shape_A[0] = 4; shape_A[1] = 3;
  shape_B[0] = 3; shape_B[1] = 2;
  shape_C[0] = shape_A[0]; shape_C[1] = shape_B[1];
  matmul(1,A,B,0,C,shape_A,shape_B,false);

  /* C = */
  /* 273     455 */
  /* 243     235 */
  /* 244     205 */
  /* 102     160 */
  cout << "a*b=C-matric from blas: \n";
  for(int i=0;i<shape_C[0];i++){
    for(int j=0;j<shape_C[1];j++){
      cout << C[i*shape_C[1]+j] << "\t";
    }
    cout << endl;
  }

  /* assert( C[0] == 273.0 );  assert( C[1] == 455.0 ); */
  /* assert( C[2] == 243.0 );  assert( C[3] == 235.0 ); */
  /* assert( C[4] == 244.0 );  assert( C[5] == 205.0 ); */
  /* assert( C[6] == 102.0 );  assert( C[7] == 160.0 ); */

  delete A;delete B; delete C;
  delete shape_A; delete shape_B; delete shape_C;
}


void matmul(double alpha, double* A,double* B, double beta,double* C,
            int* shape_A,int* shape_B ,bool transposeA){
  /* Wrapper for blas3, dgemm: Matrix/matrix multiplication. */
  /* Perform C := alpha*op( A )*op( B ) + beta*C, */

  double shape_C[2];
  shape_C[0] = shape_A[0]; shape_C[1] = shape_B[1];

  /* see /usr/include/cblas.h */
  enum CBLAS_TRANSPOSE transA;
  int lda;
  if (transposeA == true){
    transA = CblasTrans;
    lda = shape_A[0];
  }
  else{
    transA = CblasNoTrans;
    lda = shape_A[1];
  }

  /* http://www.netlib.org/blas/dgemm.f */
  /* In C arrays are stored as row-major. Fortran & Matlab use column major.
     dgemm needs the stride(LDA/LDB/LDC) of the array. Thus for RowMajor, the
     imput should be the number of columns */
  /*  Matrix A = */
  /* [1 2 3] */
  /* [4 5 6] */
  /* Row-major stores values as {1,2,3,4,5,6} */
  /* Stride here is 3 */

  /* Col-major stores values as {1, 4, 2, 5, 3, 6} */
  /* Stride here is 2 */

  cblas_dgemm(CblasRowMajor,transA , CblasNoTrans,
              shape_C[0], shape_C[1], shape_A[1], /* ==  m,n,k, */
              alpha,
              A,lda, /* Transpose: shape_A[0] */
              B,shape_B[1],
              beta,
              C,shape_C[1]);

  /* DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */
  /* Perform C := alpha*op( A )*op( B ) + beta*C, */

  /* M specifies  the number of rows of the matrix op( A ) and C */
  /* N specifies the number of columns of the matrix op( B ) and C */
  /* K specifies the number of columns of the matrix op( A ) and the number of
     rows of the matrix op( B ) */

}



void construct_Ke(MATPROPclass *mat,FEMclass *mesh,double *Ke){
  /* numerical integration of B'*C*B */


  int lengthB = mesh->sizeB[0]*mesh->sizeB[1];
  int lengthC = mesh->sizeC[0]*mesh->sizeC[1];
  double *B = new double[lengthB];
  double *C = new double[lengthC];
  double *helpMat = new double[mesh->sizeB[1]*mesh->sizeC[1]];
  double detjac = 1.0/8.0*mesh->dlx*mesh->dly*mesh->dlz;
  double factor;
  int shapeA[2] = {0};
  /* placeholder for gll points */
  double x,y,z;

  /* NOT necessary to  do memset because maple sets all values of the array.
     Also these being zero.*/
  /* memset(B,0,sizeof(double)*lengthB); */
  /* memset(C,0,sizeof(double)*lengthC); */
  /* memset(helpMat,0,sizeof(double)*mesh->sizeB[1]*mesh->sizeC[1]); */
  /* cout << "sizeB: " << lengthB << endl; */
  /* get C */
  construct_CD(mat,C);

  /* cout << "C-mat: " << endl; */
  /* for(int i=0;i<6;i++){ */
  /*   for(int j=0;j<6;j++){ */
  /*     cout << C[i*6+j] << "\t"; */
  /*   } */
  /*   cout << endl; */
  /* } */

  int tmp = mesh->nen*3;
  for(int i=0;i<mesh->nInt;i++){
    for(int j=0;j<mesh->nInt;j++){
      for(int k=0;k<mesh->nInt;k++){
        /* get B for the current gll points */
        x = mesh->intPoints[i]; y = mesh->intPoints[j]; z = mesh->intPoints[k];
        switch (mesh->ngll){
        case(3):
          construct_B3(x,y,z,mesh,B); break;
        case(4):
          construct_B4(x,y,z,mesh,B); break;
        case(5):
          construct_B5(x,y,z,mesh,B); break;
        case(6):
          construct_B6(x,y,z,mesh,B); break;
        case(7):
          construct_B7(x,y,z,mesh,B); break;
        }

        /* cout << "B-mat: " << endl; */
        /* for(int ii=0;ii<6;ii++){ */
        /*   for(int jj=0;jj<6;jj++){ */
        /*     cout << B[ii*tmp+jj] << "\t"; */
        /*   } */
        /*   cout << endl; */
        /* } */
        /* return; */
        factor = mesh->w[i]*mesh->w[j]*mesh->w[k]*detjac;
        /* helpMat = Transpose(B)*C */
        shapeA[0] = mesh->sizeB[1]; shapeA[1] = mesh->sizeB[0];
        matmul(1,B,C,0,helpMat,shapeA,mesh->sizeC,true);

        /* cout << "Help-mat: " << endl; */
        /* for(int ii=0;ii<10;ii++){ */
        /*   for(int jj=0;jj<6;jj++){ */
        /*     cout << helpMat[ii*6+jj] << "\t"; */
        /*   } */
        /*   cout << endl; */
        /* } */
        /* return; */
        
        /* Ke = Ke + helpMat*B*factor */
        shapeA[0] = mesh->sizeB[1]; shapeA[1] = mesh->sizeC[1];
        matmul(factor,helpMat,B,1,Ke,shapeA,mesh->sizeB,false);

        /* cout << "factor " << factor << "\t detjac: " << detjac << endl; */
        /* cout << "K-mat: " << endl << endl; */
        /* for(int ii=0;ii<10;ii++){ */
        /*   for(int jj=0;jj<6;jj++){ */
        /*     cout << Ke[ii*tmp+jj] << "\t"; */
        /*   } */
        /*   cout << endl; */
        /* } */
        /* return; */
      }
    }
  }
  /* cout << "K-mat: " << endl; */
  /* for(int ii=0;ii<10;ii++){ */
  /*   for(int jj=0;jj<8;jj++){ */
  /*     cout << Ke[ii*tmp+jj] << "\t"; */
  /*   } */
  /*   cout << endl; */
  /* } */

  delete[] B; delete[] C; delete[] helpMat;
}


void construct_Me(MATPROPclass *mat,FEMclass *mesh,double *M){
  /* Get diagonal mass vector */

  switch (mesh->ngll){
  case(3):
    construct_Me3(mat,mesh,M); break;
  case(4):
    construct_Me4(mat,mesh,M); break;
  case(5):
    construct_Me5(mat,mesh,M); break;
  case(6):
    construct_Me6(mat,mesh,M); break;
  case(7):
    construct_Me7(mat,mesh,M); break;
  }
}
