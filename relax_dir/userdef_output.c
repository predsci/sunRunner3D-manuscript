#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
{
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
  int i,j,k;
  double ***T;
  double ***p, ***rho;
  double mu = 0.6;

  double *x1 = grid->x[IDIR];

  T   = GetUserVar("T");


  rho = d->Vc[RHO]; /* pointer shortcut to density */
  p   = d->Vc[PRS]; /* pointer shortcut to pressure */
  
  DOM_LOOP(k,j,i){
    T[k][j][i] = (KELVIN*mu)*p[k][j][i]/rho[k][j][i];

  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





