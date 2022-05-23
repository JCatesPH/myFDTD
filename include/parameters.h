#ifndef _DOMAIN
#define _DOMAIN

#include "physicalConstants.h"
#define NUM_THREADS 6

// --- Domain size ---
#define NZ 10200      // Number of x grid points
#define NOUT 80     // Number of time steps where output is stored
#define DOUT 1000   // Spacing between time steps saved
#define LZ 200        // Num cells in PML
#define BZ 400     // Num of cells in "buffer" (SF but not PML)
#define DZ 5.0e-9     // x grid size
#define NF 1200        // Number of frequencies

// --- PML Parameters ---
#define SX_P 3          // Polynomial order
#define SX_R 1e-10      // Desired reflectivity coefficient
#define SX_M 4.0        // Max value of s_x


/* - Source parameters - */
//  Currently point source at (ic,jc)
static const int isrc     = LZ + BZ;    // Index of source location
static const double I0    = 1.0e12;  //50e16;       // Mean pump intensity [W / m^2]
static const double t0    = 0.75e-13;      // Time of pulse peak [s]
static const double taup  = 20e-15;    // Pulse duration (FWHM) for fund. pulse, second harmonic has half duration [s]
static const double lamb0 = 1.9e-6; // Vacuum wavelength of central frequency
static const double omeg0 = 2.0 * PI * C0 / lamb0;
static const double R     = 0.0;         // Ratio of fund. and 2nd harmonic
static const double phi   = 0.0;         // Phase offset of 2nd harmonic

// --- Material Parameters ---
static const int    idie1 = 3000 + LZ + BZ;    // Index where dielectric starts
static const int    idie2 = 7000 + LZ + BZ;    // Index where dielectric starts
static const double eps_inf = 2.34;     // High-frequency permittivity
static const double eps_static = 1.75;     // Low-frequency permittivity
static const double tau_Debye = 2.65e-15;     // Debye relaxation time
static const double sigma0 = 2e3; // https://link.springer.com/article/10.1007/s11468-018-0811-6
static const double n2 = 2.3e-20;      // Nonlinear index [m^2 / W]

static const double num_atoms = 2.1e26;
static const double mpi_sigmak = 1e-23; //1.0e-224; // [m^(28) W^(-14) s^(-1)] MPI cross section
static const double mpi_k = 3.0; // Photons for MPI
static const double nu_e = 1.0 / 3e-15; // Mean collision frequency for electrons
#endif
