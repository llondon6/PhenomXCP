
#ifndef _LALSIM_IMR_PHENOMX_PNR_DEVIATIONS_H
#define _LALSIM_IMR_PHENOMX_PNR_DEVIATIONS_H
/*
 * Copyright (C) 2021 The University of Amsterdam
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */


/**
 * \author Lionel London
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


// 
#include <math.h>


/* Header file for IMRPhenomXCP's tuned parameters */

// Header formatting for MU1
double IMRPhenomXCP_MU1( double theta, double eta, double a1 );
// Header formatting for MU3
double IMRPhenomXCP_MU3( double theta, double eta, double a1 );
// Header formatting for MU4
double IMRPhenomXCP_MU4( double theta, double eta, double a1 );
// Header formatting for NU4
double IMRPhenomXCP_NU4( double theta, double eta, double a1 );
// Header formatting for NU5
double IMRPhenomXCP_NU5( double theta, double eta, double a1 );
// Header formatting for NU6
double IMRPhenomXCP_NU6( double theta, double eta, double a1 );
// Header formatting for ZETA1
double IMRPhenomXCP_ZETA1( double theta, double eta, double a1 );
// Header formatting for ZETA2
double IMRPhenomXCP_ZETA2( double theta, double eta, double a1 );

#ifdef __cplusplus
}
#endif

#endif
