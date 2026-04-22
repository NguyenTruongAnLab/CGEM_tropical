#ifndef DEFINE_H
#define DEFINE_H

/* Static definitions */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define MAXM 1000       // Maximum grid points
#define MAX_TIMESERIES 50000  // Maximum time series length for forcing/BC data
#define G   9.81        // Gravitational acceleration [m/s²]
#define MAXITS 50       // Max iteration steps for tidal convergence
#define EPS   0.00001

typedef enum {
	TRIB_INJECTION_FISCHER = 0,   // Fischer (1979) mass source: S_point = Q_trib × C_trib / V_cell
	TRIB_INJECTION_RUTHERFORD = 1  // Rutherford (1994) lateral mixing
} TributaryMassInjectionMode;

/* Common includes */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>    // For clock_t and time functions


/* Global parameters */
extern long MAXT;    // Total simulation time [seconds]
extern long WARMUP;  // Warm-up period [seconds]
extern int DELTI;   // Time step [seconds]
extern int TS;      // Save every TS time step

/* Grid and geometry parameters */
extern int DELXI;   // Spatial step size [meters]
extern int EL;      // Estuarine length [meters]
extern int M;       // Max even grid points
extern int M1;      // M-1: max odd grid points
extern int M2;      // M-2: last even grid point
extern int M3;      // last odd grid point
extern double AMPL; // Tidal amplitude [m]

/* Hydrodynamic constants */
#define TIDAL_PERIOD     44712.0  /* M2 tidal period [s] */
#define MIN_DISPERSION     1.0    /* Minimum dispersion coefficient [m²/s] */

#endif // DEFINE_H
