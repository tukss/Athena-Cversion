#include "copyright.h"
/*============================================================================*/
/*! \file convbox.c
 *  \brief Problem generator for the convbox problem.
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
#error : The problem generator convbox.c does not work for mhd.
#endif /* MHD */

/*----------------------------------------------------------------------------*/
/* problem:   */

static Real interpol(Real x, int n, const Real *r, const Real *v) {
  int i;

  if (x < r[0])
    ath_error("below lower boundary for interpolation\n");

  for (i = 0; i < n; i++) {
    if (r[i] > x)
      break;
  }
  i--;
  if (i == n - 1)
    ath_error("above upper boundary for interpolation\n");

  return (x - r[i]) / (r[i+1] - r[i]) * (v[i+1] - v[i]) + v[i];
}

static Real grav_pot2(const Real x1, const Real x2, const Real x3)
{
  return 1e3 * x2;
}

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,is,ie,j,js,je,ks,nx1,nx2;
  Real d0,p0,v0,x1,x2,x3;
  int len;
  Real *r, *rho, *p, *T;

  FILE *infile;

  infile = fopen("input.dat", "r");
  if (NULL == infile)
    ath_error("input.dat missing\n");

  if (1 != fscanf(infile, "%d", &len))
    ath_error("malformed input.dat\n");

  r = calloc_1d_array(len, sizeof *r);
  if (NULL == r)
    ath_error("allocation error\n");
  rho = calloc_1d_array(len, sizeof *rho);
  if (NULL == rho)
    ath_error("allocation error\n");
  p = calloc_1d_array(len, sizeof *p);
  if (NULL == p)
    ath_error("allocation error\n");
  T = calloc_1d_array(len, sizeof *T);
  if (NULL == T)
    ath_error("allocation error\n");

  for (i = 0; i < len; i++) {
    if (4 != fscanf(infile, "%lf %lf %lf %lf", r+i, rho+i, p+i, T+i))
      ath_error("missing values in input.dat\n");
  }

  fclose(infile);

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  if (pGrid->Nx[2] > 1) {
    ath_error("[convbox]: This problem can only be run in 2D\n");
  }
  if (pGrid->Nx[1] == 1) {
    ath_error("[convbox]: This problem can only be run in 2D\n");
  }

  StaticGravPot = grav_pot2;

  v0 = 1e-9;

/* Initialize conservative fields */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
/* Calculate the cell center positions */
      cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);

      d0 = interpol(x2, len, r, rho);
      p0 = interpol(x2, len, r, p);
      pGrid->U[ks][j][i].d = d0;
      pGrid->U[ks][j][i].M1 = d0*v0*drand48();
      pGrid->U[ks][j][i].M2 = d0*v0*drand48();
      pGrid->U[ks][j][i].M3 = 0.0;
#ifndef ISOTHERMAL
    pGrid->U[ks][j][i].E = p0/Gamma_1
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
