#include "copyright.h"
/*============================================================================*/
/*! \file hotbubble.c
 *  \brief Problem generator for the hotbubble problem.
 *
 * REFERENCE: */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef MHD
#error : The problem generator hotbubble.c does not work for mhd.
#endif /* MHD */

/*----------------------------------------------------------------------------*/
/* problem:   */

static int len;
static Real *r_in, *rho_in, *p_in, *pot_in;

static Real interpol(Real x, int n, const Real *r, const Real *v) {
  int i;

  if (x < r[0])
    ath_error("below lower boundary for interpolation\n");
  else if (x > r[n-1])
    ath_error("above upper boundary for interpolation\n");

  i = (int) (n * (x - r[0]) / (r[n-1] - r[0]));
  if (i == n) i--;

  return (x - r[i]) / (r[i+1] - r[i]) * (v[i+1] - v[i]) + v[i];
}

static Real grav_pot2(const Real x1, const Real x2, const Real x3)
{
  return interpol(x2, len, r_in, pot_in);
}

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,is,ie,j,js,je,ks,nx1,nx2;
  Real rho_hse,rho,theta,dtheta,dtheta0,temp,pres,pres0,x,y,z,r,x0,y0,r0;
  const Real gasR = 8.314472e7;

  FILE *infile;

  infile = fopen("input2.dat", "r");
  if (NULL == infile)
    ath_error("input.dat missing\n");

  if (1 != fscanf(infile, "%d", &len))
    ath_error("malformed input.dat\n");

  r_in = calloc_1d_array(len, sizeof *r_in);
  if (NULL == r_in)
    ath_error("allocation error\n");
  rho_in = calloc_1d_array(len, sizeof *rho_in);
  if (NULL == rho_in)
    ath_error("allocation error\n");
  p_in = calloc_1d_array(len, sizeof *p_in);
  if (NULL == p_in)
    ath_error("allocation error\n");
  pot_in = calloc_1d_array(len, sizeof *pot_in);
  if (NULL == pot_in)
    ath_error("allocation error\n");

  for (i = 0; i < len; i++) {
    if (4 != fscanf(infile, "%lf %lf %lf %lf", r_in+i, rho_in+i, p_in+i, pot_in+i))
      ath_error("missing values in input.dat\n");
  }

  fclose(infile);

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  if (pGrid->Nx[2] > 1) {
    ath_error("[hotbubble]: This problem can only be run in 2D\n");
  }
  if (pGrid->Nx[1] == 1) {
    ath_error("[hotbubble]: This problem can only be run in 2D\n");
  }

  StaticGravPot = grav_pot2;

  x0 = 5.0e5;
  y0 = 2.75e5;
  r0 = 1.25e5;
  dtheta0 = 6.6e-4;

  // reference pressure
  pres0 = interpol(0.0, len, r_in, p_in);

/* Initialize conservative fields */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
/* Calculate the cell center positions */
      cc_pos(pGrid,i,j,ks,&x,&y,&z);

      r = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))/r0;

      if (r < 1.0) {
        dtheta = dtheta0*pow(cos(0.5 * M_PI * r), 2);
      }
      else {
        dtheta = 0.0;
      }

      rho_hse = interpol(y, len, r_in, rho_in);
      pres = interpol(y, len, r_in, p_in);

      temp = pres/(gasR*rho_hse);
      theta = temp*pow((pres0/pres), Gamma);
      theta = theta + dtheta;
      temp = theta*pow((pres/pres0), Gamma);
      rho = pres/(gasR*temp);

      pGrid->U[ks][j][i].d = rho;
      pGrid->U[ks][j][i].M1 = 0.0;
      pGrid->U[ks][j][i].M2 = 0.0;
      pGrid->U[ks][j][i].M3 = 0.0;
#ifndef ISOTHERMAL
      pGrid->U[ks][j][i].E = pres/Gamma_1
          + 0.5*(SQR(pGrid->U[ks][j][i].M1) + SQR(pGrid->U[ks][j][i].M2)
                + SQR(pGrid->U[ks][j][i].M3))/pGrid->U[ks][j][i].d;
      if (pGrid->U[ks][j][i].E <= 0.0)
        ath_error("non-positive energy\n");
#endif
    }
  }

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}
