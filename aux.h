// -*- coding: utf-8 -*-
#ifndef AUX_H
#define AUX_H

#include "fedata.h"

void matmul_test();
/* void matmul(double* A,double* B,double* C, */
/*             int* shape_A,int* shape_B,bool transposeA); */
void matmul(double alpha, double* A,double* B, double beta,double* C,
            int* shape_A,int* shape_B,bool transposeA);


void construct_Ke(MATPROPclass *mat,FEMclass *mesh,double *Ke);
void construct_B(double x, double y, double z,FEMclass *mesh, double *BB);
void construct_CD(MATPROPclass *mat, double *CCD);
void construct_Me(MATPROPclass *mat,FEMclass *mesh,double *M);
#endif
