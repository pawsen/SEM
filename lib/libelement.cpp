
#include "libelement.h"
#include <math.h>

/* char path[30]; */
/* path = "../matrix/"; */

void construct_CD(MATPROPclass *mat, double *CCD){
  /* Get elastic constitutive matrix'(C) from file */

  double E=mat->e,nu=mat->nu;

#include "../matrix/CD.c"

}

void construct_B3(double x, double y, double z,FEMclass *mesh, double *BB){
  /* Get strain/displacement ../matrix(B) file */
  /* integration points(x,y,z) are given by gll points */

  /* lengths are element lengths */
  double lx=mesh->dlx, ly=mesh->dly, lz=mesh->dlz;

#include "../matrix/B3.c"
}

void construct_B4(double x, double y, double z,FEMclass *mesh, double *BB){
  double lx=mesh->dlx, ly=mesh->dly, lz=mesh->dlz;
#include "../matrix/B4.c"
}

void construct_B5(double x, double y, double z,FEMclass *mesh, double *BB){
  double lx=mesh->dlx, ly=mesh->dly, lz=mesh->dlz;
#include "../matrix/B5.c"
}

void construct_B6(double x, double y, double z,FEMclass *mesh, double *BB){
  double lx=mesh->dlx, ly=mesh->dly, lz=mesh->dlz;
#include "../matrix/B6.c"
}

void construct_B7(double x, double y, double z,FEMclass *mesh, double *BB){
  double lx=mesh->dlx, ly=mesh->dly, lz=mesh->dlz;
#include "../matrix/B7.c"
}



void construct_Me3(MATPROPclass *mat,FEMclass *mesh,double *M){
  /* Get diagonal mass vector */

  /* lengths are element lengths */
  double lx=mesh->dlx, ly=mesh->dly, lz=mesh->dlz;
  double rho=mat->rho;

#include "../matrix/Me3.c"
}

void construct_Me4(MATPROPclass *mat,FEMclass *mesh,double *M){
  double lx=mesh->dlx, ly=mesh->dly, lz=mesh->dlz;
  double rho=mat->rho;
#include "../matrix/Me4.c"
}
void construct_Me5(MATPROPclass *mat,FEMclass *mesh,double *M){
  double lx=mesh->dlx, ly=mesh->dly, lz=mesh->dlz;
  double rho=mat->rho;
#include "../matrix/Me5.c"
}
 void construct_Me6(MATPROPclass *mat,FEMclass *mesh,double *M){
  double lx=mesh->dlx, ly=mesh->dly, lz=mesh->dlz;
  double rho=mat->rho;
#include "../matrix/Me6.c"
}
void construct_Me7(MATPROPclass *mat,FEMclass *mesh,double *M){
  double lx=mesh->dlx, ly=mesh->dly, lz=mesh->dlz;
  double rho=mat->rho;
#include "../matrix/Me7.c"
}
