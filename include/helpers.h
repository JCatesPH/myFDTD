#ifndef _HELPERS
#define _HELPERS

#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <memory>
#include <iostream>
//#include <format>
#include <fstream>
#include <string>
#include <chrono>
#include <cmath>
#include <complex.h>
#include <omp.h>

#include "parameters.h"

// --- Constants and simple macros
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
typedef std::chrono::steady_clock Clock;

// Define global variables
extern double gi2[NZ], gi3[NZ];                     // PML parameters (see Sullivan, Ch. 3)
extern double fi1[NZ], fi2[NZ], fi3[NZ];
extern double Sx[NZ], Ix[NZ];
extern double rho[NZ], Jx[NZ];
extern double Dx[NZ], Hy[NZ];   // Field arrays
extern double Ex[NZ];                           // Actual field in z-direction   
extern double chi3;
extern double Ew_re[NF][NZ], Ew_im[NF][NZ];     // Real and imaginary frequency domain output.
extern double Fsrc_re[NF], Fsrc_im[NF]; // Real and imaginary frequency domain source.
extern double eRefl_re[NF], eRefl_im[NF];
extern double omeg[NF];
extern double amp_in[NF];

// --- Function declarations ---
double src(double t);
void pmldef();      // Defines parameters for PML regions
void savegrid();    // Save XY-grid for plotting later
void initArrays(double dt);
double nfdtdsteps(int N, double T, double dt);    // Save XY-grid for plotting later
void writeSimParameters();

template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    std::unique_ptr<char[]> buf( new char[ size ] ); 
    snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

#endif
