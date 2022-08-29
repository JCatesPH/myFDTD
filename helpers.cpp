#include "include/helpers.h"

using namespace std;

/* Source function */
double src(double t) {
    return sqrt(2*I0/(C0*EPS0)) * (sqrt(1.0-R)*exp(-2.0*log(2.0)*pow((t-t0)/taup, 2)) * cos(omeg0*(t-t0)) 
        + sqrt(R) * exp(-8.0*log(2.0)*pow((t-t0)/taup, 2)) * cos(2.0*omeg0*(t-t0) + phi) );
}

double atomicDensityProfile(int index) {
	double z0 = (idie1 + idie2) / 2.0 * DZ;
    double gamma2 = pow(0.75e-6, 2);
	double n_max = 2.0e26;

    double z = index * DZ;
	//return n_max * (gamma2 / (pow(z - z0, 2) + gamma2));
    return n_max;
}

void pmldef() {
    double xn,xxn,xnum;
    /* ---- Calculate the PML parameters ---- */
    #pragma omp parallel for
    for (int i=0; i < NZ; i++){
        gi2[i] = 1.0;
        gi3[i] = 1.0;
        fi1[i] = 0.0;
        fi2[i] = 1.0;
        fi3[i] = 1.0;
    }

    for (int i=0; i <= LZ; i++){
        xnum = LZ - i;
        xxn = xnum / LZ;
        xn = 0.33*pow(xxn,SX_P);
        gi2[i] = 1.0 / (1.0 + xn);
        gi2[NZ-1-i] = 1.0 / (1.0 + xn);
        gi3[i] = (1.0 - xn) / (1.0 + xn);
        gi3[NZ-1-i] = (1.0 - xn) / (1.0 + xn);

        xxn = (xnum - 0.5) / LZ;
        xn = 0.33*pow(xxn, SX_P);
        fi1[i] = xn;
        fi1[NZ-2-i] = xn;
        fi2[i] = 1.0 / (1.0 + xn);
        fi2[NZ-2-i] = 1.0 / (1.0 + xn);
        fi3[i] = (1.0 - xn) / (1.0 + xn);
        fi3[NZ-2-i] = (1.0 - xn) / (1.0 + xn);
    }

}


void savegrid() {
    FILE *fp;
    fp = fopen("data/X.csv", "w");
    for (int i=0; i<NZ-1; i++){
        fprintf(fp, "%12.9f\n", DZ*i);
    }
    fprintf(fp, "%12.9f", DZ*(NZ-1));
    fclose(fp);

}

void initArrays(double dt) {
    ofstream freqfile;
    char valstr[20];

    #pragma omp parallel for
    for (int i=0; i<NZ; i++){
        Dx[i] = 0.0;
        Hy[i] = 0.0;
        rho[i] = 0.0;
        Jx[i] = 0.0;
        Sx[i] = 0.0;
        Ix[i] = 0.0;
        Ex1[i] = 0.0;
        sigx[i] = 0.0;
        //chi[i] = pow(eps_relDB, 2) - 1.0;
        chi[i] = 0.0;
    }

    //freqfile.open("data/omeg.csv");
    for (int k=0; k < NF; k++) {
        Fsrc_re[k] = 0.0;
        Fsrc_im[k] = 0.0;
        omeg[k] = 2*PI/((NOUT*DOUT-1)*dt) * k;
        eRefl_re[k] = 0.0;
        eRefl_im[k] = 0.0;
        //sprintf(valstr, "%12.9f", omeg[k]);
        //freqfile << valstr << endl;
    }
    //freqfile.close();
}

double nfdtdsteps(int N, double T, double dt) {
    double eSquared, eCubed;
    double gamma_Debye = (eps_static - eps_inf) * dt / tau_Debye;
    double kappa = eps_inf + gamma_Debye + sigma0 * dt / EPS0;
    double eta_k = pow(ETA0/2, mpi_k);
    for (int n=1; n<=N; n++){
        T = T + dt;
        /* - Update flux density */
        #pragma omp parallel for
        for (int i=1; i<NZ; i++){
            //Dx[i] = gi3[i] * Dx[i] + 0.5 * gi2[i] * (Hy[i-1] - Hy[i]);
            Dx[i] = gi3[i] * Dx[i] + 0.5 * gi2[i] * (Hy[i-1] - Hy[i]) - C0 * dt * Jx[i];
            //Dx[i] = Dx[i] + 0.5 * (Hy[i-1] - Hy[i]);
        }

        /* - Source - */
        Dx[isrc] = Dx[isrc] + src(T); // soft source
        //Dx[isrc] = src(T); // hard source

        /* - Calculate Ez - */
        #pragma omp parallel for
        for (int i=0; i < NZ-1; i++){
            if (i > idie1 && i < idie2) {
                eSquared = pow(ETA0 * Ex[i], 2);
                eCubed = eSquared * Ex[i];

                Ex[i] = (Dx[i] - Ix[i] - exp(-dt / tau_Debye) * Sx[i] + 2.0 * chi3 * eCubed) / (kappa + 3.0 * chi3 * eSquared); 
            }
            else {
                Ex[i] = (Dx[i] - Ix[i] - exp(-dt / tau_Debye) * Sx[i]) / kappa; // original
            }
            Sx[i] = gamma_Debye * Ex[i] + exp(-dt / tau_Debye) * Sx[i];
            Ix[i] = Ix[i] + sigma0 * dt / EPS0 * Ex[i];
        }
        /* - Calculate Fourier transform - */
        #pragma omp parallel for
        for (int k=0; k < NF; k++) {
            for (int i=1; i<NZ; i++){
                Ew_re[k][i] += cos(omeg[k]*T) * Ex[i];
                Ew_im[k][i] += sin(omeg[k]*T) * Ex[i];
            }
            if (T < 2*t0) {
                Fsrc_re[k] += cos(omeg[k]*T) * Ex[isrc];
                Fsrc_im[k] += sin(omeg[k]*T) * Ex[isrc];
            }
            if (T > 2.5*t0) {
                eRefl_re[k] += cos(omeg[k]*T) * Ex[isrc];
                eRefl_im[k] += sin(omeg[k]*T) * Ex[isrc];
            }
        }
        
        /* - Set Ez edges to 0, as part of the PML - */
        Ex[0] = 0.0;
        Ex[NZ-1] = 0.0;

        /* - Calculate Hx - */
        #pragma omp parallel for
        for (int i=0; i<NZ; i++){
            Hy[i] = fi3[i] * Hy[i] + 0.5 * fi2[i] * (Ex[i] - Ex[i+1]);
            //Hy[i] = Hy[i] + 0.5 * (Ex[i] - Ex[i+1]);
        }

        double num_atoms;
        #pragma omp parallel for
        for (int i=0; i<NZ; i++){
            num_atoms = atomicDensityProfile(i);
            /* - Calculate carrier density - */
            if (i > idie1 && i < idie2) {
                rho[i] = rho[i] + dt * mpi_sigmak * eta_k * pow(Ex[i], 2*mpi_k) * (num_atoms - rho[i]);
            }
            else {
                rho[i] = 0;
            }   

            if (rho[i] > num_atoms) {
                rho[i] = num_atoms;
            }

            /* - Calculate Jx - */
            Jx[i] = (1.0 - dt * nu_e) * Jx[i] +  dt * Q2OM * ETA0 * rho[i] * 0.5 * (Ex[i] + Ex1[i]);

            /* - Update field of previous time step - */
            Ex1[i] = Ex[i];
        }

    }

    return T;
}


double nfdtdsteps_DebyeKerr(int N, double T, double dt) {
    double eSquared;
    double tau_Kerr = 1e-15;

    for (int n=1; n<=N; n++){
        T = T + dt;
        /* - Update flux density */
        #pragma omp parallel for
        for (int i=1; i<NZ; i++){
            Dx[i] = gi3[i] * Dx[i] + 0.5 * gi2[i] * (Hy[i-1] - Hy[i]);
            //Dx[i] = gi3[i] * Dx[i] + 0.5 * gi2[i] * (Hy[i-1] - Hy[i]) - C0 * dt * Jx[i];
            //Dx[i] = Dx[i] + 0.5 * (Hy[i-1] - Hy[i]);
        }

        /* - Source - */
        Dx[isrc] = Dx[isrc] + src(T); // soft source
        //Dx[isrc] = src(T); // hard source

        for (int j = 0; j < 5; j++) {

            /* - Update value of susceptibility - */
            #pragma omp parallel for
            for (int i=0; i<NZ; i++){
                if (i > idie1 && i < idie2) {
                    eSquared = pow(ETA0 * Ex[i], 2);
                    //chi[i] = chi[i] + (dt / tau_Kerr) * (3.0 / 4.0 * chi3 * eSquared - chi[i]);
                    chi[i] = (tau_Kerr / dt * chi[i] + 3.0 / 4.0 * chi3 * eSquared) / (1.0 + tau_Kerr / dt);
                }
            }

            /* - Calculate Ez - */
            #pragma omp parallel for
            for (int i=0; i < NZ-1; i++){
                Ex[i] = Dx[i] / (1.0 + chi[i]);  
            }

        }
        
        /* - Calculate Fourier transform - */
        #pragma omp parallel for
        for (int k=0; k < NF; k++) {
            for (int i=1; i<NZ; i++){
                Ew_re[k][i] += cos(omeg[k]*T) * Ex[i];
                Ew_im[k][i] += sin(omeg[k]*T) * Ex[i];
            }
            if (T < 3*t0) {
                Fsrc_re[k] += cos(omeg[k]*T) * Ex[isrc];
                Fsrc_im[k] += sin(omeg[k]*T) * Ex[isrc];
            }
            if (T > 3.5*t0) {
                eRefl_re[k] += cos(omeg[k]*T) * Ex[isrc];
                eRefl_im[k] += sin(omeg[k]*T) * Ex[isrc];
            }
        }
        
        /* - Set Ez edges to 0, as part of the PML - */
        Ex[0] = 0.0;
        Ex[NZ-1] = 0.0;

        /* - Calculate Hx - */
        #pragma omp parallel for
        for (int i=0; i<NZ; i++){
            Hy[i] = fi3[i] * Hy[i] + 0.5 * fi2[i] * (Ex[i] - Ex[i+1]);
            //Hy[i] = Hy[i] + 0.5 * (Ex[i] - Ex[i+1]);
        }

        
    }

    return T;
}

double nfdtdsteps_JCP(int N, double T, double dt) {
    double eSquared;
    double tau_Kerr = 1e-15;

    for (int n=1; n<=N; n++){
        T = T + dt;
        /* - Update flux density */
        #pragma omp parallel for
        for (int i=1; i<NZ; i++){
            Dx[i] = gi3[i] * Dx[i] + 0.5 * gi2[i] * (Hy[i-1] - Hy[i]);
            //Dx[i] = gi3[i] * Dx[i] + 0.5 * gi2[i] * (Hy[i-1] - Hy[i]) - C0 * dt * Jx[i];
            //Dx[i] = Dx[i] + 0.5 * (Hy[i-1] - Hy[i]);
        }

        /* - Source - */
        Dx[isrc] = Dx[isrc] + src(T); // soft source
        //Dx[isrc] = src(T); // hard source

        for (int j = 0; j < 5; j++) {

            /* - Update value of susceptibility - */
            #pragma omp parallel for
            for (int i=0; i<NZ; i++){
                if (i > idie1 && i < idie2) {
                    eSquared = pow(ETA0 * Ex[i], 2);
                    chi[i] = chi[i] + (dt / tau_Kerr) * (3.0 / 4.0 * chi3 * eSquared - chi[i]);
                }
            }

            /* - Calculate Ez - */
            #pragma omp parallel for
            for (int i=0; i < NZ-1; i++){
                Ex[i] = Dx[i] / (1.0 + chi[i]);  
            }

        }
        
        /* - Calculate Fourier transform - */
        #pragma omp parallel for
        for (int k=0; k < NF; k++) {
            for (int i=1; i<NZ; i++){
                Ew_re[k][i] += cos(omeg[k]*T) * Ex[i];
                Ew_im[k][i] += sin(omeg[k]*T) * Ex[i];
            }
            if (T < 3*t0) {
                Fsrc_re[k] += cos(omeg[k]*T) * Ex[isrc];
                Fsrc_im[k] += sin(omeg[k]*T) * Ex[isrc];
            }
            if (T > 3.5*t0) {
                eRefl_re[k] += cos(omeg[k]*T) * Ex[isrc];
                eRefl_im[k] += sin(omeg[k]*T) * Ex[isrc];
            }
        }
        
        /* - Set Ez edges to 0, as part of the PML - */
        Ex[0] = 0.0;
        Ex[NZ-1] = 0.0;

        /* - Calculate Hx - */
        #pragma omp parallel for
        for (int i=0; i<NZ; i++){
            Hy[i] = fi3[i] * Hy[i] + 0.5 * fi2[i] * (Ex[i] - Ex[i+1]);
            //Hy[i] = Hy[i] + 0.5 * (Ex[i] - Ex[i+1]);
        }

        
    }

    return T;
}

std::string SIM_DATA_OUTPUT("data/");

void writeSimParameters() {
    ofstream simParamFile;
    simParamFile.open(SIM_DATA_OUTPUT + "SimParameters.dat");

    simParamFile << "NZ    \t" << NZ << "  \t/* Number of spacial grid points */" << endl;
    simParamFile << "NOUT    \t" << NOUT << "  \t/* */" << endl;
    simParamFile << "DOUT    \t" << DOUT << "  \t/* */" << endl;
    simParamFile << "LZ    \t" << LZ << "  \t/* */" << endl;
    simParamFile << "BZ    \t" << BZ << "  \t/* */" << endl;
    simParamFile << "NF    \t" << NF << "  \t/* */" << endl;

    simParamFile << "iSL    \t" << isrc << "  \t/* */" << endl;
    simParamFile << "iZ1    \t" << idie1 << "  \t/* */" << endl;
    simParamFile << "iZ2    \t" << idie2 << "  \t/* */" << endl;

    simParamFile << "DZ    \t" << DZ << "  \t/* */" << endl;
    simParamFile << "DT    \t" << DZ/6e8 << "  \t/* */" << endl;

    simParamFile << "PML_Pz    \t" << SX_P << "  \t/* */" << endl;
    simParamFile << "PML_Rz    \t" << SX_R << "  \t/* */" << endl;
    simParamFile << "PML_Mz    \t" << SX_M << "  \t/* */" << endl;

    simParamFile << "I0    \t" << I0 << "  \t/* */" << endl;
    simParamFile << "t0    \t" << t0 << "  \t/* */" << endl;
    simParamFile << "taup    \t" << taup << "  \t/* */" << endl;
    simParamFile << "lamb0    \t" << lamb0 << "  \t/* */" << endl;
    simParamFile << "omeg0    \t" << omeg0 << "  \t/* */" << endl;
    simParamFile << "Rsh    \t" << R << "  \t/* */" << endl;
    simParamFile << "phi    \t" << phi << "  \t/* */" << endl;

    simParamFile << "eps_inf    \t" << eps_inf << "  \t/* */" << endl;
    simParamFile << "eps_s    \t" << eps_static << "  \t/* */" << endl;
    simParamFile << "tau_Debye    \t" << tau_Debye << "  \t/* */" << endl;
    simParamFile << "sigma0    \t" << sigma0 << "  \t/* */" << endl;
    simParamFile << "n2    \t" << n2 << "  \t/* */" << endl;
    simParamFile << "mpi_sigmak    \t" << mpi_sigmak << "  \t/* */" << endl;
    simParamFile << "mpi_k    \t" << mpi_k << "  \t/* */" << endl;
    simParamFile << "nu_e    \t" << nu_e << "  \t/* */" << endl;

    simParamFile.close();
}


