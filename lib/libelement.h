
#ifndef LIBELEMENT_H
#define LIBELEMENT_H

#include "../fedata.h"

void construct_CD(MATPROPclass *mat, double *CCD);
void construct_B3(double x, double y, double z,FEMclass *mesh, double *BB);
void construct_B4(double x, double y, double z,FEMclass *mesh, double *BB);
void construct_B5(double x, double y, double z,FEMclass *mesh, double *BB);
void construct_B6(double x, double y, double z,FEMclass *mesh, double *BB);
void construct_B7(double x, double y, double z,FEMclass *mesh, double *BB);

void construct_Me3(MATPROPclass *mat,FEMclass *mesh,double *M);
void construct_Me4(MATPROPclass *mat,FEMclass *mesh,double *M);
void construct_Me5(MATPROPclass *mat,FEMclass *mesh,double *M);
void construct_Me6(MATPROPclass *mat,FEMclass *mesh,double *M);
void construct_Me7(MATPROPclass *mat,FEMclass *mesh,double *M);
#endif
