#ifndef GLOBALS_H
#define GLOBALS_H  
/*==============================================================================
 * FILE: globals.h
 *
 * PURPOSE: Contains following global variables:
 *     CourNo       - The Courant, Friedrichs, & Lewy (CFL) Number
 *     Iso_csound   - Isothermal sound speed
 *     Iso_csound2  - Iso_csound^2
 *     Gamma       - Gamma C_p/C_v
 *     Gamma_1, Gamma_2 - Gamma-1, and Gamma-2
 *   The first occurence in this file is included in main.c and defines the
 *   variables.  The second is included everywhere else. 
 *============================================================================*/

#ifdef MAIN_C

Real CourNo;
#ifdef ISOTHERMAL
Real Iso_csound;
Real Iso_csound2;
#else
Real Gamma;
Real Gamma_1, Gamma_2;
#endif
StaticGravAcc_t x1GravAcc = NULL, x2GravAcc = NULL, x3GravAcc = NULL;

/*----------------------------------------------------------------------------*/
/* definitions included everywhere except main.c  */

#else /* MAIN_C */

extern Real CourNo;
#ifdef ISOTHERMAL
extern Real Iso_csound;
extern Real Iso_csound2;
#else
extern Real Gamma;
extern Real Gamma_1, Gamma_2;
#endif
extern StaticGravAcc_t x1GravAcc, x2GravAcc, x3GravAcc;

#endif /* MAIN_C */
#endif /* GLOBALS_H */