/* From the Sullivan FDTD book */
#include "include/parameters.h"
#include "include/helpers.h"

using namespace std;

// Define global variables
double gi2[NZ], gi3[NZ];                     // PML parameters (see Sullivan, Ch. 3)
double fi1[NZ], fi2[NZ], fi3[NZ];            // PML parameters (see Sullivan, Ch. 3) 
double Sx[NZ], Ix[NZ];
double rho[NZ], Jx[NZ];                   // Carrier and current density arrays
double Dx[NZ], Hy[NZ];           // Field arrays
double Ex[NZ], Ex1[NZ];                           // Actual field in z-direction   
double chi3;
double Ew_re[NF][NZ], Ew_im[NF][NZ]; // Real and imaginary frequency domain output.
double Fsrc_re[NF], Fsrc_im[NF];     // Real and imaginary frequency domain source.
double eRefl_re[NF], eRefl_im[NF];
double omeg[NF];                     // Frequencies for DFT.
double amp_in[NF];


int main()
{
    auto t1 = Clock::now();

    omp_set_num_threads(NUM_THREADS);
	cout << "Num threads set to  = " << NUM_THREADS << endl << "  Test: ";

	#pragma omp parallel 
	{
		#pragma omp critical
		{
			cout << omp_get_thread_num() << " ";
		}
	}
	cout << endl << endl;

    /* ---- Define variables ---- */
    int n,i,j; // indices
    string fname;

    /* - Integration parameters - */
    int nsteps[NOUT]; // Number of time steps
    double dt = DZ / 6e8; // Time step (dx/2*c) 
    double T = 0.0; // Initial time

    for (n=0; n<NOUT; n++){
        nsteps[n] = n*DOUT;
    }

    chi3 = 4 * EPS0 * C0 * n2 / 3;

    /* ---- Print information ---- */
    cout << "FDTD parameters: " << endl;
    cout << "  Number of time steps : ";
    for (n=0; n<NOUT; n++){
        printf("%4d ", nsteps[n]);
    }
    cout << endl << "-----------------------------" << endl;
    cout << " Grid parameters:" << endl;
    cout << "   Num. grid     : " << NZ << endl;
    cout << "   Cells in PML  : " << LZ << endl;
    cout << "   Buffer region : " << BZ << endl;
    cout << "   Grid spacing  : " << DZ << endl;
    cout << "   Time Step     : " << dt << endl;
    cout << "   Total time    : " << NOUT*DOUT*dt << endl;
    cout << "-----------------------------" << endl;
    cout << " PML parameters:" << endl;
    cout << "   Polynomial    : " << SX_P << endl;
    cout << "   Reflectivity  : " << SX_R << endl;
    cout << "   Smax          : " << SX_M << endl;
    cout << "-----------------------------" << endl;
    cout << " Beam parameters: (two-color Gaussian pulse)" << endl;
    cout << "   Source index  : " << isrc << endl;
    cout << "   Intensity     : " << I0 << endl;
    cout << "   Pulse peak    : " << t0 << endl;
    cout << "   FWHM dur.     : " << taup << endl;
    cout << "   Wavelength    : " << lamb0 << endl;
    cout << "   Frequency     : " << omeg0 << endl;
    cout << "   2nd Ratio     : " << R << endl;
    cout << "   Phase offset  : " << phi << endl;
    /* cout << "-----------------------------" << endl;
    cout << " Material parameters:" << endl;
    cout << "   Rel. Perm.    : " << epsz << endl;
    cout << "   Conductivity  : " << sigma << endl;
    cout << "   Nonlinear ind.: " << n2 << endl;
    cout << "   Chi3          : " << chi3 << endl; */


    writeSimParameters();

    /* - Save x-y grid - */
    savegrid();
    /* - Set PML variables - */
    pmldef();

    /* ---- Initialize arrays ---- */
    initArrays(dt);

    char valstr[50];
    ofstream outfile;
    outfile.open("data/omeg.csv");
    for (int k=0; k < NF; k++) {
        sprintf(valstr, "%10.7e", omeg[k]);
        outfile << valstr << endl;
    }
    outfile.close();

    /* ---- Main FDTD Loop ---- */
    printf("=============================\n");
    printf("Starting FDTD simulation. \n");
    printf("=============================\n");
    printf("  n  |  nsteps(n)  |  T  \n");
    printf("-----------------------------\n");

    ofstream eField_file, rho_file, reflRe_outfile, reflIm_outfile;
    char re_outname[50], im_outname[50];
    eField_file.open("data/Ex.csv");
    rho_file.open("data/rho.csv");
    /* sourceRe_outfile.open();
    sourceIm_outfile.close(); */
    
    for (n=0; n<NOUT; n++){

        printf("%3d  |  %5d  | %8.5e\n", n, nsteps[n], T);
        T = nfdtdsteps(DOUT, T, dt);
        
        for (i=0; i<NZ-1; i++){
            sprintf(valstr, "%12.9f, ", Ex[i]);
            eField_file << valstr;

            sprintf(valstr, "%12.9e, ", rho[i]);
            rho_file << valstr;
        }
        sprintf(valstr, "%12.9f", Ex[NZ-1]);
        eField_file << valstr << endl;

        sprintf(valstr, "%12.9e", rho[i]);
        rho_file << valstr << endl;

        /* -- Print the frequency information to files*/
        /* sprintf(re_outname, "data/Ew_re_%d.csv", n*DOUT);
        sprintf(im_outname, "data/Ew_re_%d.csv", n*DOUT);
        reflRe_outfile.open(re_outname);
        reflIm_outfile.open(im_outname);
        for (int i=0; i<NZ; i++){
            for (int k=0; k < NF-1; k++) {
                //sprintf(valstr, "%12.9e, ", Ew_re[k][i] / amp_in[k]);
                sprintf(valstr, "%12.9e, ", Ew_re[k][i]);
                reflRe_outfile << valstr;
                //reflRe_outfile << string_format("%12.9e, ", Ew_re[k][i]);
                //sprintf(valstr, "%12.9e, ", Ew_im[k][i] / amp_in[k]);
                sprintf(valstr, "%12.9e, ", Ew_im[k][i]);
                reflIm_outfile << valstr;
                //reflIm_outfile << string_format("%12.9e, ", Ew_im[k][i]);
            }
            //sprintf(valstr, "%12.9e", Ew_re[NF-1][i] / amp_in[NF-1]);
            sprintf(valstr, "%12.9e\n", Ew_re[NF-1][i]);
            reflRe_outfile << valstr;
            //sprintf(valstr, "%12.9e", Ew_im[NF-1][i] / amp_in[NF-1]);
            sprintf(valstr, "%12.9e\n", Ew_im[NF-1][i]);
            reflIm_outfile << valstr;
        }
        reflRe_outfile.close();
        reflIm_outfile.close();  */
    }
    eField_file.close();    
    rho_file.close();
    
    printf("=============================\n");
    /* --- End of main loop --- */
    
    reflRe_outfile.open("data/Ew_re.csv");
    reflIm_outfile.open("data/Ew_im.csv");
    for (int i=0; i<NZ; i++){
        for (int k=0; k < NF-1; k++) {
            sprintf(valstr, "%12.9e, ", Ew_re[k][i]);
            reflRe_outfile << valstr;
            sprintf(valstr, "%12.9e, ", Ew_im[k][i]);
            reflIm_outfile << valstr;
        }
        sprintf(valstr, "%12.9e\n", Ew_re[NF-1][i]);
        reflRe_outfile << valstr;
        sprintf(valstr, "%12.9e\n", Ew_im[NF-1][i]);
        reflIm_outfile << valstr;
    }
    reflRe_outfile.close();
    reflIm_outfile.close();

    reflRe_outfile.open("data/eRefl_re.csv");
    reflIm_outfile.open("data/eRefl_im.csv");
    for (int k=0; k < NF-1; k++) {
        sprintf(valstr, "%12.9e, ", eRefl_re[k]);
        reflRe_outfile << valstr;
        sprintf(valstr, "%12.9e, ", eRefl_im[k]);
        reflIm_outfile << valstr;
    }
    sprintf(valstr, "%12.9e\n", eRefl_re[NF-1]);
    reflRe_outfile << valstr;
    sprintf(valstr, "%12.9e\n", eRefl_im[NF-1]);
    reflIm_outfile << valstr;
    reflRe_outfile.close();
    reflIm_outfile.close();

    reflRe_outfile.open("data/eSrc_re.csv");
    reflIm_outfile.open("data/eSrc_im.csv");
    for (int k=0; k < NF-1; k++) {
        sprintf(valstr, "%12.9e, ", Fsrc_re[k]);
        reflRe_outfile << valstr;
        sprintf(valstr, "%12.9e, ", Fsrc_im[k]);
        reflIm_outfile << valstr;
    }
    sprintf(valstr, "%12.9e\n", Fsrc_re[NF-1]);
    reflRe_outfile << valstr;
    sprintf(valstr, "%12.9e\n", Fsrc_im[NF-1]);
    reflIm_outfile << valstr;
    reflRe_outfile.close();
    reflIm_outfile.close();

    auto t2 = Clock::now();
    printf("Simulation finished successfully. \n");
    cout << "  Time to calculate: " 
        <<  (double) std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() * 1e-6
        << " [seconds]" << std::endl;
    return 0;
}
