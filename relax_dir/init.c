/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sepy 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  double  g_omegarot = g_inputParam[OMEGA_ROT];

  v[TRC] = 0.0;
#if ROTATING_FRAME == YES
  g_OmegaZ  = g_omegarot;
#endif    
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
  int i,j,k,id;
  double ***rho;
  double alpha, beta;
  double T;
  double two = 2.0;
  double mu = 0.6;
  double g_gamma = g_inputParam[GAMMA];  
  double g_r0 = g_inputParam[R0];
	  
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
   
  beta =  two;

  alpha = two * g_gamma - two;
  
  id = InputDataOpen ("2D/vr_r0.flt","2D/vr_r0_grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i) {
    d->Vc[VX1][k][j][i] = InputDataInterpolate (id,g_r0,x2[j],x3[k]);   
    d->Vc[VX2][k][j][i] = 0.0;
    d->Vc[VX3][k][j][i] = 0.0;
    
  }
  InputDataClose(id);

  id = InputDataOpen ("2D/rho_r0.flt","2D/rho_r0_grid.out"," ", 0, CENTER);
  TOT_LOOP(k,j,i) {
    d->Vc[RHO][k][j][i] = InputDataInterpolate (id,g_r0,x2[j],x3[k]) * pow(g_r0/x1[i], beta);
  }
  InputDataClose(id);
  
  /* need to convert temperature to pressure see top of page 46 in PLUTO manual */
      
  rho   = d -> Vc[RHO]; /* pointer shortcut to density */
      
  id = InputDataOpen("2D/t_r0.flt","2D/t_r0_grid.out"," ",0, CENTER);
  TOT_LOOP(k,j,i) {
    T = InputDataInterpolate (id,g_r0,x2[j],x3[k]) * pow(g_r0/x1[i], alpha);
    d->Vc[PRS][k][j][i] = T * rho[k][j][i] /(KELVIN * mu);
  }
  InputDataClose(id);
       
#if PHYSICS == MHD || PHYSICS == RMHD
      
  id = InputDataOpen("2D/br_r0.flt","2D/br_r0_grid.out"," ",0, CENTER);
  TOT_LOOP(k,j,i) {
    d->Vc[BX1][k][j][i] = InputDataInterpolate (id,g_r0,x2[j],x3[k]) * pow(g_r0/x1[i], beta);
    d->Vc[BX2][k][j][i] = 0.0;
    d->Vc[BX3][k][j][i] = 0.0;
  }
  InputDataClose(id);
          
#endif
  
  
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * Follow the example from the Mach_reflection example.....
 *
 *********************************************************************** */
{
  int i,j,k,id, nv;

  double T;
  double two = 2.0;
  double mu = 0.6;
  double  g_gamma = g_inputParam[GAMMA];  
  double  g_omegarot = g_inputParam[OMEGA_ROT];
  double g_r0 = g_inputParam[R0];
  

  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

  double twopi = 2.0 * CONST_PI;
  
  static int first_call=1;
  
  static Data_Arr Vc0;

  double ***rho_bcr0;

  /* read in the BC files */

  if (first_call) {

    id = InputDataOpen ("2D/vr_r0.flt","2D/vr_r0_grid.out"," ", 0, CENTER);
    BOX_LOOP(box, k,j,i) {
    
      d->Vc[VX1][k][j][i] = InputDataInterpolate (id,g_r0,x2[j],x3[k]);   
      d->Vc[VX2][k][j][i] = 0.0;
      d->Vc[VX3][k][j][i] = 0.0;
    
    }
    InputDataClose(id);
   
    id = InputDataOpen("2D/rho_r0.flt","2D/rho_r0_grid.out"," ",0, CENTER);
    BOX_LOOP(box, k,j,i) {

      d->Vc[RHO][k][j][i] = InputDataInterpolate (id,g_r0,x2[j],x3[k]);
    }
    InputDataClose(id);
  
    rho_bcr0 = d -> Vc[RHO];

    /* need to convert temperature to pressure see top of page 46 in PLUTO manual */
      
    id = InputDataOpen("2D/t_r0.flt","2D/t_r0_grid.out"," ",0, CENTER);
    BOX_LOOP(box, k,j,i) {
      T = InputDataInterpolate (id,g_r0,x2[j],x3[k]);
      d->Vc[PRS][k][j][i] = T * rho_bcr0[k][j][i] /(KELVIN * mu);
    }
    InputDataClose(id);
    

#if PHYSICS == MHD || PHYSICS == RMHD
      
    id = InputDataOpen("2D/br_r0.flt","2D/br_r0_grid.out"," ",0, CENTER);
    BOX_LOOP(box, k,j,i) {
      d->Vc[BX1][k][j][i] = InputDataInterpolate (id,g_r0,x2[j],x3[k]);
      d->Vc[BX2][k][j][i] = 0.0;
      d->Vc[BX3][k][j][i] = 0.0;
    }
    InputDataClose(id);

#endif
    Vc0 = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
    TOT_LOOP(k,j,i) {
       NVAR_LOOP(nv) {
       Vc0[nv][k][j][i] = d -> Vc[nv][k][j][i];
      }
    }
    first_call = 0;
      
  } else {

    //TOT_LOOP(k,j,i) {
	BOX_LOOP(box,k,j,i){	  
      NVAR_LOOP(nv) {
          d -> Vc[nv][k][j][i] = Vc0[nv][k][j][i];
      }
    }    
  }   

}


#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
