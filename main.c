/*
 *  PoroViscoElastic wave propagation in 2D isotropic media, for one attenuation
 *  mechanism, on a staggered grid
 *
 *  Bernard Giroux
 *  INRS-ETE
 *
 *
 *  Code specs:  - language: ANSI C99
 *               - compiled with intel compiler 11.1 on a mac running OSX 10.6
 *               - external library: fftw ( http://www.fftw.org )
 *  
 *  
 *   Reference papers
 *
 
 @ARTICLE{carcione99b,
 author = {Jos\'{e} M. Carcione and Hans B. Helle},
 title = {Numerical Solution of the Poroviscoelastic Wave Equation on a Staggered Mesh},
 journal = {Journal of Computational Physics},
 year = {1999},
 volume = {154},
 pages = {520 - 527},
 number = {2},
 doi = {10.1006/jcph.1999.6321}
 }
 
 @ARTICLE{carcione96d,
 author = {Jos\'{e} M. Carcione},
 title = {Wave propagation in anisotropic, saturated porous media: Plane-wave
 theory and numerical simulation},
 journal = {Journal of the Acoustical Society of America},
 year = {1996},
 volume = {99},
 pages = {2655-2666},
 number = {5},
 doi = {10.1121/1.414809}
 }
 
 @ARTICLE{carcione95b,
 author = {Jos\'e M. Carcione and Gerardo Quiroga-Goode},
 title = {Some aspects of the physics and numerical modeling of {B}iot compressional waves},
 journal = {Journal of Computational Acoustics},
 year = {1995},
 volume = {3},
 pages = {261--280},
 number = {4},
 doi = {10.1142/S0218396X95000136}
 } 
 *
 *
 */

/*
 
 ------------------------------------------------------------------------------

 Usage
 
 pve_iso [options] -m model.in -s source.in -o output_description.in

   options are:
     -t value    with "value" the length of the time window in s (1.2 by default)
     -d value    with "value" the time step in ms (0.1 by default)
     -h Print this message
     -v verbose mode (type twice for increased verbosity)
 
 
 
 ------------------------------------------------------------------------------

 Format of model file
 
 --- beginning of file ---
 2 2            <-  nx and nz, number of nodes in x & z
 1.0 1.0        <-  dx and dz, grid step size in x & z, in meters
 0.0 0.0        <-  origin of the grid, in meters
 K_m(1,1) K_s(1,1) K_f(1,1) phi(1,1) mu(1,1) rho_s(1,1) rho_f(1,1) T(1,1) eta(1,1) kappa(1,1) Q(1,1)
 K_m(1,2) K_s(1,2) K_f(1,2) phi(1,2) mu(1,2) rho_s(1,2) rho_f(1,2) T(1,2) eta(1,2) kappa(1,2) Q(1,2)
 K_m(2,1) K_s(2,1) K_f(2,1) phi(2,1) mu(2,1) rho_s(2,1) rho_f(2,1) T(2,1) eta(2,1) kappa(2,1) Q(2,1)
 K_m(2,2) K_s(2,2) K_f(2,2) phi(2,2) mu(2,2) rho_s(2,2) rho_f(2,2) T(2,2) eta(2,2) kappa(2,2) Q(2,2)
 --- end of file ---
 
 Variable  Description                                 Units
 
 K_m       bulk modulus of drained matrix            [ GPa ]
 K_s       bulk modulus of the solid                 [ GPa ]
 K_f       bulk modulus of the fluid                 [ GPa ]
 phi       porosity                                    [ - ]
 mu        shear modulus of the matrix               [ GPa ]
 rho_s     solid density                          [ kg/m^3 ]
 rho_f     fluid density                          [ kg/m^3 ]
 T         tortuosity                                  [ - ]
 eta       fluid viscosity                            [ cP ]
 kappa     permeability                               [ mD ]
 Q         seismic quality factor                      [ - ]
 
 note: 1 cP = 1e-3 Pa.s ; 1 mD = 9.869233e-16 m^2
 
 
 
 ------------------------------------------------------------------------------
 
 Format of source file
 
 --- beginning of file ---
 2              <- number of source points
 
 Sx             <- Source component of 1st source point
 10.0           <- Strength in MPa at 1st source point
 50.0           <- Frequency in Hz at 1st source point
 50.0 45.0      <- x z coordinates at 1st source point in meters
 
 Sf             <- Source component of 2nd source point
 1.0            <- Strength in MPa at 2nd source point
 50.0           <- Frequency in Hz at 2nd source point
 50.0 45.0      <- x z coordinates at 2nd source point in meters
 --- end of file ---
 
 Components can be Sx, Sz, Sxz, Sf or Bulk
 
 
 
 ------------------------------------------------------------------------------
 
 Format of file defining the outputs
 
 --- beginning of file ---
 2              <- number of records
 test1          <- basename of all output files
 
 Vx             <- particle velocity component of 1st record
 Trace          <- type of 1st record
 10.0 20.0      <- X Z coordinates at 1st record point in meters
 0.0            <- start time of 1st recording, in seconds
 2.0            <- time sampling of 1st record, in miliseconds
 
 Qx             <- particle velocity component of 2nd record
 Snapshot       <- type of 2nd record
 0.0            <- start time of snapshot(s) (2nd recording), in seconds
 2.0            <- time sampling of snapshot(s) 2nd record, in miliseconds.
                     Give 0.0 or less for just one snapshot at time given above.
 --- end of file ---

 All output files are ascii files.  Files names are
 
 basename_NO.trc            for traces, where NO is the record number
 basename_NO_TIME.ssh       for snapshots, where NO is record number and TIME is the record time
 
 Trace files are two columns, first for time and second for velocity component.
   Velocity components in traces are _not_ interpolated at exact X and Z coordinates,
   but rather taken at the closest grid point.
 Snapshot files are three columns, respectively for X, Z and velocity component.
 
 Units of particle velocity are mm/s, time in s and distances in m.
 
 */

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "io_utils.h"
#include "propagate.h"
#include "structs.h"

int verbose=0;

int main (int argc, char *argv[]) {

    struct inputParams params;
    struct computationVariables c;
    struct materialProperties mp;
    struct fftw_data fd;
	struct sourceParams src;
    struct outputParams out;
    
    // -------------------------------------------------------------------------
    //
    // variables
    //
    // -------------------------------------------------------------------------

    const double pi = 4.0*atan(1.0);
    const double L = 1.0;  // number of attenuation mechanisms
    const char *typeSrc[] = { "Sx", "Sz", "Sxz", "Sf", "Bulk" };
    const char *component[] = { "Vx", "Vz", "qx", "qz" };

    //
    // time variables
    //
    size_t nsteps=0;
    
    //
    // Grid properties
    //
    struct grid g;
    
    //
    // Source parameters
    //
    double f0=0.0;         // nominal frequency
    

    //
    // Motion variables
    //
    double *tau_xx;    // stress in x on plane normal to x
    double *tau_zz;    // stress in z on plane normal to z
    double *tau_xz;    // stress in x on plane normal to z
    double *p;         // fluid pressure
    double *v_x;       // solid particle velocity along x
    double *v_z;       // solid particle velocity along z
    double *q_x;       // fluid particle velocity (relative to solid) along x
    double *q_z;       // fluid particle velocity (relative to solid) along z

    double *vs_x;      // solid particle velocity along x
    double *vs_z;      // solid particle velocity along z
    double *qs_x;      // fluid particle velocity (relative to solid) along x
    double *qs_z;      // fluid particle velocity (relative to solid) along z
    
    double *tau_xxtmp;
    double *tau_zztmp;
    double *tau_xztmp;
    double *ptmp;

    double *e;         // memory variable (one att. mechanism), 
    double *m;         // term in Darcy's law
    double *rho;       // composite density                           [ kg/m^3 ]
    
    double *coeff_v_i; // coefficients for analytical solution of the stiff system
    double *coeff_v_j;    
    double *coeff_q_i;
    double *coeff_q_j;
    
    double *W;
    double *Ws;
    double *Wtmp;
    double *delta;
    
    FILE *fid;
    
    set_defaults(&g, &params);
    process_args(argc, argv, &params);
    c.dt = params.dt;
    
    if ( verbose ) {
        fprintf(stdout, "\n *** pve_iso - PoroViscoElastic wave propagation in 2D isotropic media ***\n\n");
        if ( verbose > 1 ) {
            fprintf(stdout, "  Model file: \t%s\n", params.modelfile);
            fprintf(stdout, "  Source file:\t%s\n", params.sourcefile);
            fprintf(stdout, "  Output file:\t%s\n", params.outputfile);
            fprintf(stdout, "  Time window:\t%lg (s)\n", params.duration);
            fprintf(stdout, "  Time step:  \t%lg (ms)\n", params.dt * 1.e3);
        }
        fprintf(stdout, "Reading size of model ... ");
    }
    read_grid_params(params.modelfile, &g);
    
    if ( verbose ) {
        fprintf(stdout, "done.\n");
        if ( verbose > 1 ) {
            fprintf(stdout, "  Model parameters: \t\tX\t\tZ\n");
            fprintf(stdout, "        grid points \t\t%zd\t\t%zd\n", g.nx, g.nz);
            fprintf(stdout, "        grid spacing\t\t%lg\t\t%lg\n", g.dx, g.dz);
            fprintf(stdout, "        grid origin \t\t%lg\t\t%lg\n", g.x0, g.z0);
            fprintf(stdout, "        grid padding\t\t%zd\t\t%zd\n", g.npx, g.npz);
        }
        fprintf(stdout, "Allocating memory ... ");
        fflush(stdout);
    }
    
    g.nx2 = g.nx + 2*g.npx;
    g.nz2 = g.nz + 2*g.npz;

    //
    //  Memory allocation
    //
    
    size_t nnodes = g.nx2 * g.nz2;
    
    if ( NULL == ( mp.K_m   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.K_s   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.K_f   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.phi   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.mu    = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.rho_s = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.rho_f = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.T     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.eta   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.kappa = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.Q     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }

    if ( NULL == ( c.epsilon = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.E       = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.M       = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.alpha   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.varphi  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.tau_s   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.rho_i   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.rho_j   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.rho_f_i = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.rho_f_j = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.nk_i    = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.nk_j    = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.mu_ij   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.m_i     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.m_j     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }

    if ( NULL == ( m         = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( rho       = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
    if ( NULL == ( coeff_v_i = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( coeff_v_j = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( coeff_q_i = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( coeff_q_j = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
    c.mu = &(mp.mu[0]);
    
    //
    // Read in model
    // 
    
    if ( verbose ) {
        fprintf(stdout, "done.\nReading properties of materials ... ");
        fflush(stdout);
    }
    read_model(params.modelfile, &g, &mp);
    
    //
    // Read in source params
    //
    if ( verbose ) {
        fprintf(stdout, "done.\nReading source parameters ... ");
        fflush(stdout);
    }
    read_source(params.sourcefile, &g, &src);
	f0 = src.s[0].f;
	compute_src_fct(&src, params.dt);
	
    //
    // Read in output params
    //
    if ( verbose ) {
        fprintf(stdout, "done.\n");
        if ( verbose > 1 ) {
            fprintf(stdout, "  Source has %zd component(s)\n", src.nsrc);
            for ( size_t ns=0; ns<src.nsrc; ++ns ) {
                fprintf(stdout, "    %zd - type:       \t%s\n", ns+1, typeSrc[src.s[ns].type]);
                fprintf(stdout, "    %zd - frequency:  \t%lg Hz\n", ns+1, src.s[ns].f);
                fprintf(stdout, "    %zd - strength:   \t%lg MPa\n", ns+1, src.s[ns].A*1.e3);
                fprintf(stdout, "    %zd - coordinates:\t(%lg, %lg)\n", ns+1, src.s[ns].x, src.s[ns].z);
            }
        }
        fprintf(stdout, "Reading output parameters ... ");
        fflush(stdout);
    }
    read_output(params.outputfile, &g, &out);
    
    if ( verbose ) {
        fprintf(stdout, "done.\n");
        if ( verbose > 1 ) {
            fprintf(stdout, "  There will be %zd record(s) on output\n", out.nrec);
            for ( size_t nr=0; nr<out.nrec; ++nr ) {
                if ( out.r[nr].type == TRACE ) {
                    fprintf(stdout, "    %zd - Trace of %s recorded at (%lf, %lf), ",
                            nr+1, component[out.r[nr].comp], out.r[nr].x, out.r[nr].z);
                    fprintf(stdout, "starting at %lg s and sampled every %lg ms\n",
                            out.r[nr].t0, 1.e3*out.r[nr].dt);
                }
                else if ( out.r[nr].type == SNAPSHOT ) {
                    if ( out.r[nr].dt <= 0.0 )
                        fprintf(stdout, "    %zd - Snapshot of %s at time %lg s\n",
                                nr+1, component[out.r[nr].comp], out.r[nr].t0);
                    else
                        fprintf(stdout, "    %zd - Snapshot of %s, starting at %lg s and sampled every %lg ms\n",
                                nr+1, component[out.r[nr].comp], out.r[nr].t0, 1.e3*out.r[nr].dt);
                }
            }
        }
        fprintf(stdout, "Preprocessing data ... ");
        fflush(stdout);
    }
    
	//
    //  Scale input data
	//
	for ( size_t n=0; n<nnodes; ++n ) {
		mp.kappa[n] *= 9.869233e-4;    // eta/kappa will be in GPa.s/m^2
		mp.rho_s[n] *= 1.e-9;          // density are now in Tg/m^3 (Giga. kg/m^3)
		mp.rho_f[n] *= 1.e-9;
    }
//	for ( size_t n=0; n<nnodes; ++n ) {
//        mp.K_m[n] *= 1.e-9;
//        mp.K_s[n] *= 1.e-9;
//        mp.K_f[n] *= 1.e-9;
//        mp.mu[n]  *= 1.e-9;
//        mp.eta[n] *= 1.e-3;
//        mp.kappa[n] *= 9.869233e-16;
//	}
//    for (size_t ns=0; ns<src.nsrc; ++ns)
//        for (size_t ne=0; ne<src[ns].length; ++ne)
//            src[ns].fct[ne] *= 1.e-9;
	
	//
    // Compute elastic & other physical params
    //
    double tau0 = 1./(2. * pi * f0);
    for ( size_t n=0; n<nnodes; ++n ) {
        c.E[n] = mp.K_m[n] + 4./3.*mp.mu[n];                                 // eq (6)
        double D = mp.K_s[n] * (1. + mp.phi[n]*(mp.K_s[n]/mp.K_f[n] - 1. )); // eq (8)
        c.M[n] = mp.K_s[n]*mp.K_s[n] / (D - mp.K_m[n]);                      // eq (7)
        c.alpha[n] = 1. - mp.K_m[n]/mp.K_s[n];                               // eq (9)
        
        m[n] = mp.T[n]*mp.rho_f[n] / mp.phi[n];
        
        double tau_e = tau0/mp.Q[n] * (sqrt(mp.Q[n]*mp.Q[n] + 1.) + 1.);
        c.tau_s[n]   = tau0/mp.Q[n] * (sqrt(mp.Q[n]*mp.Q[n] + 1.) - 1.);
        c.varphi[n]  = tau_e/c.tau_s[n] - 1.;                                // eq (19)
        
        rho[n] = (1.-mp.phi[n])*mp.rho_s[n] + mp.phi[n]*mp.rho_f[n];
    }
    
    //
    //  Averaging of parameters at staggered nodes, at (i+1/2,j)
    //
    for ( size_t i=0; i<g.nx2-1; ++i ) {
        for ( size_t j=0; j<g.nz2; ++j ) {
            
            size_t ij = i*g.nz2+j;
            size_t iij = (i+1)*g.nz2+j;
            
            c.rho_i[ij]   = 0.5 * ( rho[ij]      + rho[iij] );
            c.rho_f_i[ij] = 0.5 * ( mp.rho_f[ij] + mp.rho_f[iij] );
            c.nk_i[ij]    = 0.5 * ( mp.eta[ij]/mp.kappa[ij] +
                                   mp.eta[iij]/mp.kappa[iij] );
            c.m_i[ij]     = 0.5 * ( m[ij]        + m[iij] );
            
            // Coefficients for analytical solution, at (i+1/2,j)
            double tmp = 1. / (c.rho_f_i[ij]*c.rho_f_i[ij] -
                               c.rho_i[ij]*c.m_i[ij]);
            double beta12 = tmp * c.rho_f_i[ij];                // eq (42) of Carcione (1996)
            double beta22 = tmp * -c.rho_i[ij];
            double lambda_s = -(c.nk_i[ij]) * beta22;           // eq (71) of Carcione (1996)
            coeff_q_i[ij] = exp(lambda_s*c.dt);                 // eq (77) of Carcione (1996)
            coeff_v_i[ij] = beta12/beta22*(coeff_q_i[ij] - 1.); // eq (76) of Carcione (1996)
        }
    }
    
    //
    //  Averaging of parameters at staggered nodes, at (i,j+1/2)
    //
    for ( size_t i=0; i<g.nx2; ++i ) {
        for ( size_t j=0; j<g.nz2-1; ++j ) {
            
            size_t ij = i*g.nz2+j;
            size_t ijj = ij+1;

            c.rho_j[ij]   = 0.5 * ( rho[ij]      + rho[ijj] );
            c.rho_f_j[ij] = 0.5 * ( mp.rho_f[ij] + mp.rho_f[ijj] );
            c.nk_j[ij]    = 0.5 * ( mp.eta[ij]/mp.kappa[ij] +
                                   mp.eta[ijj]/mp.kappa[ijj] );
            c.m_j[ij]     = 0.5 * ( m[ij]        + m[ijj] );
            
            // Coefficients for analytical solution, at (i,j+1/2)
            double tmp = 1. / (c.rho_f_j[ij]*c.rho_f_j[ij] -
                               c.rho_j[ij]*c.m_j[ij]);
            double beta12 = tmp * c.rho_f_j[ij];                // eq (42) of Carcione (1996)
            double beta22 = tmp * -c.rho_j[ij];
            double lambda_s = -(c.nk_j[ij]) * beta22;           // eq (71) of Carcione (1996)
            coeff_q_j[ij] = exp(lambda_s*c.dt);                 // eq (77) of Carcione (1996)
            coeff_v_j[ij] = beta12/beta22*(coeff_q_j[ij] - 1.); // eq (76) of Carcione (1996)
        }
    }
    
    //
    //  Averaging of parameters at staggered nodes, at (i+1/2,j+1/2)
    //
    for ( size_t i=0; i<g.nx2-1; ++i ) {
        for ( size_t j=0; j<g.nz2-1; ++j ) {
            
            size_t ij = i*g.nz2+j;
    
            c.mu_ij[ij] = 4.0 / ( 1./mp.mu[ij] +
                                  1./mp.mu[(i+1)*g.nz2+j] +
                                  1./mp.mu[i*g.nz2+j+1] +
                                  1./mp.mu[(i+1)*g.nz2+j+1] );
        }            
    }
    //  fill end of grid
    for ( size_t i=g.nx2-1, j=0; j<g.nz2; ++j ) {
        
        size_t ij = (i-1)*g.nz2+j;
        size_t iij = i*g.nz2+j;
        
        c.rho_i[iij]     = c.rho_i[ij];
        c.rho_f_i[iij]   = c.rho_f_i[ij];
        c.nk_i[iij]      = c.nk_i[ij];
        c.m_i[iij]       = c.m_i[ij];
        
        coeff_q_i[iij]   = coeff_q_i[ij];
        coeff_v_i[iij]   = coeff_v_i[ij];
        
        c.mu_ij[iij]     = c.mu_ij[ij];
    }
    for ( size_t i=0, j=g.nz2-1; i<g.nx2; ++i ) {
        
        size_t ij = i*g.nz2+j-1;
        size_t ijj = i*g.nz2+j;
        
        c.rho_j[ijj]     = c.rho_j[ij];
        c.rho_f_j[ijj]   = c.rho_f_j[ij];
        c.nk_j[ijj]      = c.nk_j[ij];
        c.m_j[ijj]       = c.m_j[ij];
        
        coeff_q_j[ijj]   = coeff_q_j[ij];
        coeff_v_j[ijj]   = coeff_v_j[ij];

        c.mu_ij[ijj]     = c.mu_ij[ij];
    }
    
	
	
    //
    // For effective computation:
    // storing   1/tau_s                                into   tau_s
    //           M*varphi/(1+varphi)                    into   varphi
    //           E-2mu                                  into   mu
    //           rho_f^2/rho - m                        into   m
    //           rho_f/rho                              into   rho_f
    //           1/rho                                  into   rho
    //
    for ( size_t n=0; n<nnodes; ++n ) {
        c.tau_s[n] = 1./c.tau_s[n];
        c.varphi[n] = c.M[n]*c.varphi[n]/(1.+c.varphi[n]);
        
        c.mu[n] = c.E[n] - 2.*c.mu[n];
        
        c.m_i[n] = c.rho_f_i[n]*c.rho_f_i[n]/c.rho_i[n] - c.m_i[n];
        c.m_j[n] = c.rho_f_j[n]*c.rho_f_j[n]/c.rho_j[n] - c.m_j[n];
        
        c.rho_f_i[n] /= c.rho_i[n];
        c.rho_f_j[n] /= c.rho_j[n];
        
        c.rho_i[n] = 1./c.rho_i[n];
        c.rho_j[n] = 1./c.rho_j[n];
    }
	
	
	free ( mp.Q );
    free ( mp.kappa );
    free ( mp.eta );
    free ( mp.T );
    free ( mp.rho_f );
    free ( mp.rho_s );
    free ( mp.phi );
    free ( mp.K_f );
    free ( mp.K_s );
    free ( mp.K_m );
	// mp.mu not freed yet because accessed by c struct component
	free ( rho );
    free ( m );

    if ( verbose ) {
        fprintf(stdout, "done.\nAllocating (large) motion variable vectors ... ");
        fflush(stdout);
    }
    
    
	if ( NULL == ( W      = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( Ws     = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( Wtmp   = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( delta  = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	
    for ( size_t n=0; n<9*nnodes; ++n )
        W[n] = Ws[n] = delta[n] = 0.0;
	
    tau_xx = &(W[0]);
    tau_zz = &(W[nnodes]);
    tau_xz = &(W[2*nnodes]);
    p      = &(W[3*nnodes]);
    v_x    = &(W[4*nnodes]);
    v_z    = &(W[5*nnodes]);
    q_x    = &(W[6*nnodes]);
    q_z    = &(W[7*nnodes]);
    e      = &(W[8*nnodes]);
    
    vs_x    = &(Ws[4*nnodes]);
    vs_z    = &(Ws[5*nnodes]);
    qs_x    = &(Ws[6*nnodes]);
    qs_z    = &(Ws[7*nnodes]);
	
    tau_xxtmp = &(Wtmp[0]);
    tau_zztmp = &(Wtmp[nnodes]);
    tau_xztmp = &(Wtmp[2*nnodes]);
    ptmp      = &(Wtmp[3*nnodes]);
	
    if ( verbose ) {
        fprintf(stdout, "done.\nInitializing fftw data (this may take some time) ... ");
        fflush(stdout);
    }
    init_fftw_data(Wtmp, &fd, &g);

        
	//
	//  Main loop
	//
    
    nsteps = params.duration / c.dt;
    if ( nsteps*c.dt < params.duration ) nsteps++;
    
    if ( verbose )
        fprintf(stdout, "done.\nStarting main loop (%zd iterations)\n", nsteps);
    for ( size_t it=0; it<nsteps; ++it ) {
        
		if ( (it+1)%100 == 0 && verbose ) printf("  Iteration %zd\n", it+1);
		
        // analytical solutions
        for ( size_t n=0; n<nnodes; ++n ) {
            // eq (76) & eq (77) of Carcione (1996)
            
            vs_x[n] = v_x[n] + coeff_v_i[n]*q_x[n];
            vs_z[n] = v_z[n] + coeff_v_j[n]*q_z[n];
            qs_x[n] = coeff_q_i[n]*q_x[n];
            qs_z[n] = coeff_q_j[n]*q_z[n];
        }
        // update remaining components of Ws
        for ( size_t n=0; n<4*nnodes; ++n )
            Ws[n] = W[n];
        for ( size_t n=8*nnodes; n<9*nnodes; ++n )
            Ws[n] = W[n];
        
        // eq (80) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n )
            Wtmp[n] = Ws[n];  // -> using Wtmp rather that Ws for first ride because fftw plans use Wtmp
        propagate(Wtmp, delta, &g, &c, &fd);

		
		
		
        // eq (81) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt/2. * delta[n];
			W[n]    = Ws[n] + c.dt/6. * delta[n];
		}
		// add source term
		for (size_t ns=0; ns<src.nsrc; ++ns) {
			src.s[ns].it = 2*it;
			if ( (src.s[ns].type == SF || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				ptmp[ src.s[ns].i ] += c.dt/2. * src.s[ns].fct[ src.s[ns].it ];
				p[ src.s[ns].i ]    += c.dt/6. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( (src.s[ns].type == SX || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				tau_xxtmp[ src.s[ns].i ] += c.dt/2. * src.s[ns].fct[ src.s[ns].it ];
				tau_xx[ src.s[ns].i ]    += c.dt/6. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( (src.s[ns].type == SZ || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				tau_zztmp[ src.s[ns].i ] += c.dt/2. * src.s[ns].fct[ src.s[ns].it ];
				tau_zz[ src.s[ns].i ]    += c.dt/6. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( src.s[ns].type == SXZ && src.s[ns].it<src.s[ns].length ) {
				tau_xztmp[ src.s[ns].i ] += c.dt/2. * src.s[ns].fct[ src.s[ns].it ];
				tau_xz[ src.s[ns].i ]    += c.dt/6. * src.s[ns].fct[ src.s[ns].it ];
			}
		}
        propagate(Wtmp, delta, &g, &c, &fd);        

        
		
        // eq (82) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt/2. * delta[n];
			W[n]   +=         c.dt/3. * delta[n];
		}
		// add source term
		for (size_t ns=0; ns<src.nsrc; ++ns) {
			src.s[ns].it = 2*it+1;
			if ( (src.s[ns].type == SF || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				ptmp[ src.s[ns].i ] += c.dt/2. * src.s[ns].fct[ src.s[ns].it ];
				p[ src.s[ns].i ]    += c.dt/3. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( (src.s[ns].type == SX || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				tau_xxtmp[ src.s[ns].i ] += c.dt/2. * src.s[ns].fct[ src.s[ns].it ];
				tau_xx[ src.s[ns].i ]    += c.dt/3. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( (src.s[ns].type == SZ || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				tau_zztmp[ src.s[ns].i ] += c.dt/2. * src.s[ns].fct[ src.s[ns].it ];
				tau_zz[ src.s[ns].i ]    += c.dt/3. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( src.s[ns].type == SXZ && src.s[ns].it<src.s[ns].length ) {
				tau_xztmp[ src.s[ns].i ] += c.dt/2. * src.s[ns].fct[ src.s[ns].it ];
				tau_xz[ src.s[ns].i ]    += c.dt/3. * src.s[ns].fct[ src.s[ns].it ];
			}
		}
        propagate(Wtmp, delta, &g, &c, &fd);
        
		
		
        // eq (83) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt    * delta[n];
			W[n]   +=         c.dt/3. * delta[n];
		}
		// add source term
		for (size_t ns=0; ns<src.nsrc; ++ns) {
			src.s[ns].it = 2*it+1;
			if ( (src.s[ns].type == SF || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				ptmp[ src.s[ns].i ] += c.dt    * src.s[ns].fct[ src.s[ns].it ];
				p[ src.s[ns].i ]    += c.dt/3. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( (src.s[ns].type == SX || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				tau_xxtmp[ src.s[ns].i ] += c.dt    * src.s[ns].fct[ src.s[ns].it ];
				tau_xx[ src.s[ns].i ]    += c.dt/3. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( (src.s[ns].type == SZ || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				tau_zztmp[ src.s[ns].i ] += c.dt    * src.s[ns].fct[ src.s[ns].it ];
				tau_zz[ src.s[ns].i ]    += c.dt/3. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( src.s[ns].type == SXZ && src.s[ns].it<src.s[ns].length ) {
				tau_xztmp[ src.s[ns].i ] += c.dt    * src.s[ns].fct[ src.s[ns].it ];
				tau_xz[ src.s[ns].i ]    += c.dt/3. * src.s[ns].fct[ src.s[ns].it ];
			}
		}
        propagate(Wtmp, delta, &g, &c, &fd);
        
        

        // eq (79) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n )
            W[n] += c.dt/6. * delta[n];
		// add source term
		for (size_t ns=0; ns<src.nsrc; ++ns) {
			src.s[ns].it = 2*it+2;
			if ( (src.s[ns].type == SF || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				p[ src.s[ns].i ]      += c.dt/6. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( (src.s[ns].type == SX || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				tau_xx[ src.s[ns].i ] += c.dt/6. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( (src.s[ns].type == SZ || src.s[ns].type == BULK) && src.s[ns].it<src.s[ns].length ) {
				tau_zz[ src.s[ns].i ] += c.dt/6. * src.s[ns].fct[ src.s[ns].it ];
			}
			if ( src.s[ns].type == SXZ && src.s[ns].it<src.s[ns].length ) {
				tau_xz[ src.s[ns].i ] += c.dt/6. * src.s[ns].fct[ src.s[ns].it ];
			}
		}
		
		
		for ( size_t nr=0; nr<out.nrec; ++nr ) {
            double t = c.dt*(it+1);
            double *data;
            switch (out.r[nr].comp) {
                case VX:
                    data = v_x;
                    break;
                case VZ:
                    data = v_z;
                    break;
                case QX:
                    data = q_x;
                    break;
                case QZ:
                    data = q_z;
                    break;
                default:
                    break;
            }
            if ( ( fabs(out.r[nr].t0-t) < c.dt ) ||
                ( ( t > out.r[nr].t0 && out.r[nr].dt > 0.0 &&
                   ( fmod((t-out.r[nr].t0),out.r[nr].dt) < 1.e-10 ) ) ) ) {
                if ( out.r[nr].type == TRACE )
                    write_trace(data, t, &out, nr);
                else if ( out.r[nr].type == SNAPSHOT )
                    write_snapshot(data, t, &g, &out, nr);
            }
        }
    }
	
	for ( size_t n=0; n<out.nrec; ++n )
        if ( out.r[n].type == TRACE )
            fclose( out.r[n].fid );

	for (size_t n=0; n<src.nsrc; ++n) free( src.s[n].fct );
	free ( src.s );
    free_fftw_data(&fd);

    free ( delta );
    free ( Wtmp );
    free ( Ws );
    free ( W );
    free ( coeff_q_j );
    free ( coeff_q_i );
    free ( coeff_v_j );
    free ( coeff_v_i );
    free ( c.m_j );
    free ( c.m_i );
    free ( c.mu_ij );
    free ( c.nk_j );
    free ( c.nk_i );
    free ( c.rho_f_j );
    free ( c.rho_f_i );
    free ( c.rho_j );
    free ( c.rho_i );
    free ( c.tau_s );
    free ( c.varphi );
    free ( c.alpha );
    free ( c.M );
    free ( c.E );
    free ( c.epsilon );
    free ( mp.mu );
    
    return 0;
}
