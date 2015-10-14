/*
 *  PoroViscoElastic wave propagation in 2D VTI anisotropic media, for one 
 *  attenuation mechanism, on a staggered grid
 *
 *  Created by Bernard Giroux on 11-01-06.
 *
 *  Code specs:  - language: ANSI C99
 *               - compiled with intel compiler 11.1 on a mac running OSX 10.6
 *               - external libraries: fftw ( http://www.fftw.org ), 
 *                       netCDF ( http://www.unidata.ucar.edu/software/netcdf/ )
 */  

/*
 * Copyright (c) 2011, Bernard Giroux
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of California, Berkeley nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE REGENTS AND CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*  
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
 
 @ARTICLE{carcione98b,
 author = {Carcione, J. M.},
 title = {Viscoelastic effective rheologies for modelling wave propagation
 in porous media},
 journal = {Geophysical Prospecting},
 year = {1998},
 volume = {46},
 pages = {249--270},
 number = {3},
 doi = {10.1046/j.1365-2478.1998.00087.x},
 publisher = {Blackwell Publishing Ltd},
 } 
 *
 *
 *
 */

/*
 
 ------------------------------------------------------------------------------
 
 Usage
 
 pve_vti [options] -p parameter_file.dat
 
 options are:
 -h Print this message
 -v verbose mode (type twice for increased verbosity)
 
 Format of parameter file:
 
 One line per parameter, first word is the value of the parameter, next is a
 keyword comprise between a # and a comma, and a can be added after the comma,
 i.e we have on each line
 
 value     # keyword, optional comment
 
 Example showing available keywords:
 --- beginning of file ---
 model.dat          # model, model file
 source.dat         # source, source file
 output.dat         # output, description of output records
 run1               # basename, common for all output files
 1.2                # time, length of time window in s (default is 1.2)
 0.1                # dt, time step in ms (default is 0.1)
 1                  # abs, apply absorbing boundary (0 or 1, 0 by default)
 20                 # nl, number of absorbing layers (20 by default)
 1                  # plstr, plot absorbing strips in snapshots
 --- end of file ---
 
 
 
 ------------------------------------------------------------------------------
 
 Format of model file
 
 --- beginning of file ---
 2 2                  <-  nx and nz, number of nodes in x & z
 1.0 1.0              <-  dx and dz, grid step size in x & z, in meters
 0.0 0.0              <-  origin of the grid, in meters
 anisotropic_vti      <-  keyword for type of model (isotropic, anisotropic_vti)
 K_m(1,1) K_s(1,1) K_f(1,1) phi(1,1) mu(1,1) rho_s(1,1) rho_f(1,1) T3(1,1) eta(1,1) kappa3(1,1) Q(1,1) f0(1,1) epsilon(1,1) delta(1,1) k3k1(1,1) T3T1(1,1)
 K_m(1,2) K_s(1,2) K_f(1,2) phi(1,2) mu(1,2) rho_s(1,2) rho_f(1,2) T3(1,2) eta(1,2) kappa3(1,2) Q(1,2) f0(1,2) epsilon(1,2) delta(1,2) k3k1(1,2) T3T1(1,2)
 K_m(2,1) K_s(2,1) K_f(2,1) phi(2,1) mu(2,1) rho_s(2,1) rho_f(2,1) T3(2,1) eta(2,1) kappa3(2,1) Q(2,1) f0(2,1) epsilon(2,1) delta(2,1) k3k1(2,1) T3T1(2,1)
 K_m(2,2) K_s(2,2) K_f(2,2) phi(2,2) mu(2,2) rho_s(2,2) rho_f(2,2) T3(2,2) eta(2,2) kappa3(2,2) Q(2,2) f0(2,2) epsilon(2,2) delta(2,2) k3k1(2,2) T3T1(2,2)
 --- end of file ---
 
 Variable  Description                                                     Units
 
 K_m       bulk modulus of drained matrix                                [ GPa ]
 K_s       bulk modulus of the solid                                     [ GPa ]
 K_f       bulk modulus of the fluid                                     [ GPa ]
 phi       porosity                                                        [ - ]
 mu        shear modulus of the matrix                                   [ GPa ]
 rho_s     solid density                                              [ kg/m^3 ]
 rho_f     fluid density                                              [ kg/m^3 ]
 T3        tortuosity (along z)                                            [ - ]
 eta       fluid viscosity                                                [ cP ]
 kappa3    permeability  (along z)                                        [ mD ]
 Q         seismic quality factor                                          [ - ]
 f0        relaxation frequency                                           [ Hz ]
 epsilon   Thomsen anisotropy parameter                                    [ - ]
 delta     Thomsen anisotropy parameter                                    [ - ]
 k3k1      permeability anisotropy ratio ( kappa_z / kappa_x )             [ - ]
 T3T1      tortuosity anisotropy ratio ( T_z / T_x )                       [ - ]
 
 note: 1 cP = 1e-3 Pa.s ; 1 mD = 9.869233e-10 m^2
 
 
 
 ------------------------------------------------------------------------------
 
 Format of source file
 
 --- beginning of file ---
 2              <- number of source points
 1              <- Number of pts in source template, can be 1 or 9 (Lin et Thylén, 2009)
 
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
 
 @ARTICLE{lin09,
 author = {Zhili Lin and Lars Thylén},
 title = {An analytical derivation of the optimum source patterns for the pseudospectral
 time-domain method},
 journal = {Journal of Computational Physics},
 year = {2009},
 volume = {228},
 pages = {7375 - 7387},
 number = {19},
 doi = {10.1016/j.jcp.2009.06.033}
 }
 
 
 ------------------------------------------------------------------------------
 
 Format of file defining the outputs
 
 --- beginning of file ---
 2              <- number of records
 
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
 
 Trace files are ascii files and snapshots are binary grid files in netCDF format (COARDS-compliant).
 Files names are
 
 basename_NO.trc            for traces, where NO is the record number
 basename_NO_TIME.nc        for snapshots, where NO is record number and TIME is the record time
 
 Trace files are two columns, first for time and second for velocity component.
 Velocity components in traces are _not_ interpolated at exact X and Z coordinates,
 but rather taken at the closest grid point.
 Snapshot grids can be viewed in matlab, with Paraview ( http://www.paraview.org/ ) or by 
 using GMT routines ( http://www.soest.hawaii.edu/gmt/ ) to create plots.
 
 Units of particle velocity are mm/s, time in s and distances in m.
 
 */

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bc.h"
#include "io_utils.h"
#include "propagate.h"
#include "structs.h"
#include "src.h"

int verbose=0;

int main (int argc, char *argv[]) {

    struct inputParams params;
    struct computationVariablesVTI c;
    struct materialPropertiesVTI mp;
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
    const char *typeSrc[] = { "Sx", "Sy", "Sz", "Sxy", "Sxz", "Syz", "Bulk", "Sf" };
    const char *component[] = { "Vx", "Vy", "Vz", "qx", "qz", "tau_xx", "tau_zz", "tau_xz", "p", "tau_xy", "tau_yz" };
	
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
    double *m1;        // term in Darcy's law
    double *m3;        // term in Darcy's law
    double *rho;       // composite density                           [ kg/m^3 ]
	
	double *coeff_v_1; // coefficients for analytical solution of the stiff system
    double *coeff_v_3;    
    double *coeff_q_1;
    double *coeff_q_3;
    
	double *W;
    double *Ws;
    double *Wtmp;
    double *delta;
	
    params.bc = &g.bc;
	set_defaults(&g, &params);
    process_args(argc, argv, &params);
    c.dt = params.dt;
	strcpy(out.basename, params.basename);

    if ( verbose ) {
        fprintf(stdout, "\n *** pve_vti - PoroViscoElastic wave propagation in 2D VTI anisotropic media ***\n\n");
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
            fprintf(stdout, "        grid padding\t\t%zd\t\t%zd\n", g.bc.np, g.bc.np);
        }
        fprintf(stdout, "Allocating memory ... ");
        fflush(stdout);
    }
    
    g.nx2 = g.nx + 2*g.bc.np;
    g.nz2 = g.nz + 2*g.bc.np;
	
	//
    //  Memory allocation
    //
    
    size_t nnodes = g.nx2 * g.nz2;
    
    if ( NULL == ( mp.K_m     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.K_s     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.K_f     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.phi     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.mu      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.rho_s   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.rho_f   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.T1      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.T3      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.eta     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.kappa1  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.kappa3  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.Q       = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.f0      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.epsilon = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.delta   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
//    if ( NULL == ( mp.gamma   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }

    if ( NULL == ( c.c11     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.c13     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.c33     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.c55     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.epsilon = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.M       = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.alpha1  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.alpha3  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.varphi  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.tau_s   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.rho_i   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.rho_j   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.rho_f_i = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.rho_f_j = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.nk1     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.nk3     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.m1      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.m3      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }

    if ( NULL == ( m1        = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( m3        = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( rho       = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
	if ( NULL == ( coeff_v_1 = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( coeff_v_3 = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( coeff_q_1 = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( coeff_q_3 = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	
	
	//
    // Read in model
    // 
    
    if ( verbose ) {
        fprintf(stdout, "done.\nReading properties of materials ... ");
        fflush(stdout);
    }
    read_modelVTI(params.modelfile, &g, &mp);
    
	if ( params.iwipe ) {
        if ( verbose ) fprintf(stdout, "done.\nApplying absorbing boundaries");
        gb2( &g, 0.92 );
    }
    else {
        if ( verbose ) fprintf(stdout, "done.\nNo absorbing boundaries");
    }
    
	//
    // Read in source params
    //
    if ( verbose ) {
        fprintf(stdout, "\nReading source parameters ... ");
        fflush(stdout);
    }
    read_source(params.sourcefile, &g, &src);
	f0 = src.s[0].f;
	compute_src_fct(&src, params.dt, 0.5);
	normalize_src_fct(&src, &g);

    //
    // Read in output params
    //
    if ( verbose ) {
        fprintf(stdout, "done.\n");
        if ( verbose > 1 ) {
            fprintf(stdout, "  Source has %zd component(s)\n", src.nsrc);
			double fac =1.;
			if ( src.nTemplate == 9 ) {
				fac = 2.;
				fprintf(stdout, "  Each component uses a %zd points template\n", src.nTemplate);
			}
            for ( size_t ns=0; ns<src.nsrc; ++ns ) {
                fprintf(stdout, "    %zd - type:       \t%s\n", ns+1, typeSrc[src.s[ns].type]);
                fprintf(stdout, "    %zd - frequency:  \t%lg Hz\n", ns+1, src.s[ns].f);
                fprintf(stdout, "    %zd - strength:   \t%lg MPa\n", ns+1, src.s[ns].A*fac);
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
            fprintf(stdout, "There will be %zd record(s) on output\n", out.nrec);
            fprintf(stdout, "  Snapshots will ");
            if ( !params.plotStrips ) fprintf(stdout, "not ");
            fprintf(stdout, "include absorbing strips.\n");
            for ( size_t nr=0; nr<out.nrec; ++nr ) {
                if ( out.r[nr].type == TRACE ) {
                    fprintf(stdout, "  %zd - Trace of %s at (%lg, %lg), ",
                            nr+1, component[out.r[nr].comp], out.r[nr].x, out.r[nr].z);
                    fprintf(stdout, "starting at %lg s & sampled at %lg ms\n",
                            out.r[nr].t0, 1.e3*out.r[nr].dt);
                }
                else if ( out.r[nr].type == SNAPSHOT ) {
                    if ( out.r[nr].dt <= 0.0 )
                        fprintf(stdout, "  %zd - Snapshot of %s at time %lg s\n",
                                nr+1, component[out.r[nr].comp], out.r[nr].t0);
                    else
                        fprintf(stdout, "  %zd - Snapshot of %s, starting at %lg s & sampled at %lg ms\n",
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
		mp.K_m[n] *= 1000.;                // now in MPa
		mp.K_s[n] *= 1000.;                // now in MPa
		mp.K_f[n] *= 1000.;                // now in MPa
		mp.mu[n]  *= 1000.;                // now in MPa
        mp.eta[n] *= 1.013249965828145e6;  // eta/kappa will be in Mega. kg/m^3/s
		mp.rho_s[n] *= 1.e-6;              // density are now in Gg/m^3 (Mega. kg/m^3)
		mp.rho_f[n] *= 1.e-6;		
    }
	
	//
    // Compute elastic tensor components & other physical params
    //
    for ( size_t n=0; n<nnodes; ++n ) {
		double tau0 = 1./(2. * pi * mp.f0[n]);
        double D = mp.K_s[n] * (1. + mp.phi[n]*(mp.K_s[n]/mp.K_f[n] - 1. )); // eq (12), carcione96d
        c.c33[n] = mp.K_m[n] + 4./3.*mp.mu[n];                               // eq (14), carcione96d
		double c55 = mp.mu[n];                                               // eq (15), carcione96d
		c.c11[n] = 2. * mp.epsilon[n] * c.c33[n] + c.c33[n];
		//double c12 = c.c11[n] - 2. * ( 2. * mp.gamma[n] * c55 + c55 );
		double c12 = c.c11[n] - 2. * c55;
		c.c13[n] = sqrt( 2. * mp.delta[n] * c.c33[n] * (c.c33[n] - c55) +
						(c.c33[n] - c55) * (c.c33[n] - c55) ) - c55;
		
		
        c.M[n] = mp.K_s[n]*mp.K_s[n] /
		     (D - ( 2.*c.c11[n]+c.c33[n] + 2.*c12 + 4.*c.c13[n])/9.);        // eq (11)
        c.alpha1[n] = 1. - 1./3. * (c.c11[n]+c12+c.c13[n])/mp.K_s[n];        // eq (9)
        c.alpha3[n] = 1. - 1./3. * (2.*c.c13[n] + c.c33[n])/mp.K_s[n];       // eq (10)
		
        m1[n] = mp.T1[n]*mp.rho_f[n] / mp.phi[n];
        m3[n] = mp.T3[n]*mp.rho_f[n] / mp.phi[n];
        
		if ( mp.Q[n] < 1.e4 ) {
			double tau_e = tau0/mp.Q[n] * (sqrt(mp.Q[n]*mp.Q[n] + 1.) + 1.);
			c.tau_s[n]   = tau0/mp.Q[n] * (sqrt(mp.Q[n]*mp.Q[n] + 1.) - 1.);
			c.varphi[n]  = tau_e/c.tau_s[n] - 1.;                                // eq (19)
		} else {
			c.tau_s[n]  = tau0;
			c.varphi[n] = 0.0;
		}
        
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
            c.nk1[ij]     = 0.5 * ( mp.eta[ij]/mp.kappa1[ij] +
                                   mp.eta[iij]/mp.kappa1[iij] );
            c.m1[ij]      = 0.5 * ( m1[ij]        + m1[iij] );
            
            // Coefficients for analytical solution, at (i+1/2,j)
            double tmp = 1. / (c.rho_f_i[ij]*c.rho_f_i[ij] -
                               c.rho_i[ij]*c.m1[ij]);
            double beta12 = tmp * c.rho_f_i[ij];                // eq (42) of Carcione (1996)
            double beta22 = tmp * -c.rho_i[ij];
            double lambda_s = -(c.nk1[ij]) * beta22;            // eq (71) of Carcione (1996)
            coeff_q_1[ij] = exp(lambda_s*c.dt);                 // eq (77) of Carcione (1996)
            coeff_v_1[ij] = beta12/beta22*(coeff_q_1[ij] - 1.); // eq (76) of Carcione (1996)
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
            c.nk3[ij]     = 0.5 * ( mp.eta[ij]/mp.kappa3[ij] +
                                   mp.eta[ijj]/mp.kappa3[ijj] );
            c.m3[ij]     = 0.5 * ( m3[ij]        + m3[ijj] );
            
            // Coefficients for analytical solution, at (i,j+1/2)
            double tmp = 1. / (c.rho_f_j[ij]*c.rho_f_j[ij] -
                               c.rho_j[ij]*c.m3[ij]);
            double beta12 = tmp * c.rho_f_j[ij];                // eq (42) of Carcione (1996)
            double beta22 = tmp * -c.rho_j[ij];
            double lambda_s = -(c.nk3[ij]) * beta22;           // eq (71) of Carcione (1996)
            coeff_q_3[ij] = exp(lambda_s*c.dt);                 // eq (77) of Carcione (1996)
            coeff_v_3[ij] = beta12/beta22*(coeff_q_3[ij] - 1.); // eq (76) of Carcione (1996)
        }
    }
    
	//
    //  Averaging of parameters at staggered nodes, at (i+1/2,j+1/2)
    //
    for ( size_t i=0; i<g.nx2-1; ++i ) {
        for ( size_t j=0; j<g.nz2-1; ++j ) {
            
            size_t ij = i*g.nz2+j;
			
            c.c55[ij] = 4.0 / ( 1./mp.mu[ij] +
							    1./mp.mu[(i+1)*g.nz2+j] +
							    1./mp.mu[i*g.nz2+j+1] +
							    1./mp.mu[(i+1)*g.nz2+j+1] );
        }            
    }
    //  fill end of grid
    for ( size_t i=g.nx2-1, j=0; j<g.nz2; ++j ) {
        
        size_t ij  = (i-1)*g.nz2+j;
        size_t iij = i*g.nz2+j;
        
        c.rho_i[iij]    = c.rho_i[ij];
        c.rho_f_i[iij]  = c.rho_f_i[ij];
        c.nk1[iij]      = c.nk1[ij];
        c.m1[iij]       = c.m1[ij];
        
        coeff_q_1[iij]  = coeff_q_1[ij];
        coeff_v_1[iij]  = coeff_v_1[ij];
        
        c.c55[iij]      = c.c55[ij];
    }
    for ( size_t i=0, j=g.nz2-1; i<g.nx2; ++i ) {
        
        size_t ij  = i*g.nz2+j-1;
        size_t ijj = i*g.nz2+j;
        
        c.rho_j[ijj]    = c.rho_j[ij];
        c.rho_f_j[ijj]  = c.rho_f_j[ij];
        c.nk3[ijj]      = c.nk3[ij];
        c.m3[ijj]       = c.m3[ij];
        
        coeff_q_3[ijj]  = coeff_q_3[ij];
        coeff_v_3[ijj]  = coeff_v_3[ij];
		
        c.c55[ijj]      = c.c55[ij];
    }
	
	//
    // For effective computation:
    // storing   1/tau_s                                into   tau_s
    //           M*varphi/(1+varphi)                    into   varphi
    //           rho_f^2/rho - m                        into   m
    //           rho_f/rho                              into   rho_f
    //           1/rho                                  into   rho
    //
    for ( size_t n=0; n<nnodes; ++n ) {
        c.tau_s[n] = 1./c.tau_s[n];
        c.varphi[n] = c.M[n]*c.varphi[n]/(1.+c.varphi[n]);
        
        c.m1[n] = c.rho_f_i[n]*c.rho_f_i[n]/c.rho_i[n] - c.m1[n];
        c.m3[n] = c.rho_f_j[n]*c.rho_f_j[n]/c.rho_j[n] - c.m3[n];
        
        c.rho_f_i[n] /= c.rho_i[n];
        c.rho_f_j[n] /= c.rho_j[n];
        
        c.rho_i[n] = 1./c.rho_i[n];
        c.rho_j[n] = 1./c.rho_j[n];
    }
	
//	free ( mp.gamma );
	free ( mp.delta );
	free ( mp.epsilon );
	free ( mp.f0 );
	free ( mp.Q );
    free ( mp.kappa3 );
    free ( mp.kappa1 );
    free ( mp.eta );
    free ( mp.T3 );
    free ( mp.T1 );
    free ( mp.rho_f );
    free ( mp.rho_s );
	free ( mp.mu );
    free ( mp.phi );
    free ( mp.K_f );
    free ( mp.K_s );
    free ( mp.K_m );
	free ( rho );
    free ( m3 );
    free ( m1 );
	
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
	//    e      = &(W[8*nnodes]);
    
    vs_x    = &(Ws[4*nnodes]);
    vs_z    = &(Ws[5*nnodes]);
    qs_x    = &(Ws[6*nnodes]);
    qs_z    = &(Ws[7*nnodes]);
	
    tau_xxtmp = &(Wtmp[0]);
    tau_zztmp = &(Wtmp[nnodes]);
    tau_xztmp = &(Wtmp[2*nnodes]);
    ptmp      = &(Wtmp[3*nnodes]);
	
    if ( verbose ) {
        fprintf(stdout, "done.\n");
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
        
		if ( params.iwipe ) {
            tbc2( &g, W );
        }
        
		double t = c.dt*(it+1);
		if ( (it+1)%50 == 0 && verbose ) { fprintf(stdout, "  Iteration %zd, t = %g s\n", it+1, t); fflush(stdout); }
		
        // analytical solutions
        for ( size_t n=0; n<nnodes; ++n ) {
            // eq (76) & eq (77) of Carcione (1996)
            
            vs_x[n] = v_x[n] + coeff_v_1[n]*q_x[n];
            vs_z[n] = v_z[n] + coeff_v_3[n]*q_z[n];
            qs_x[n] = coeff_q_1[n]*q_x[n];
            qs_z[n] = coeff_q_3[n]*q_z[n];
        }
        // update remaining components of Ws
        for ( size_t n=0; n<4*nnodes; ++n )
            Ws[n] = W[n];
        for ( size_t n=8*nnodes; n<9*nnodes; ++n )
            Ws[n] = W[n];
        
        // eq (80) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n )
            Wtmp[n] = Ws[n];  // -> using Wtmp rather that Ws for first ride because fftw plans use Wtmp
        
        propagateVTI(Wtmp, delta, &g, &c, &fd);
        
        // eq (81) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt/2. * delta[n];
			W[n]    = Ws[n] + c.dt/6. * delta[n];
		}
		// add source term
		add_src_rk4(&src, W, Wtmp, c.dt/6., c.dt/2., it, nnodes);
		
        propagateVTI(Wtmp, delta, &g, &c, &fd);        
		
        // eq (82) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt/2. * delta[n];
			W[n]   +=         c.dt/3. * delta[n];
		}
		// add source term
		add_src_rk4(&src, W, Wtmp, c.dt/3., c.dt/2., it, nnodes);
		
        propagateVTI(Wtmp, delta, &g, &c, &fd);
        
        // eq (83) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt    * delta[n];
			W[n]   +=         c.dt/3. * delta[n];
		}
		// add source term
		add_src_rk4(&src, W, Wtmp, c.dt/3., c.dt, it, nnodes);
		
        propagateVTI(Wtmp, delta, &g, &c, &fd);
        
        // eq (79) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n )
            W[n] += c.dt/6. * delta[n];
		// add source term
		add_src_rk4(&src, W, Wtmp, c.dt/6., 0.0, it, nnodes);
		
		for ( size_t nr=0; nr<out.nrec; ++nr ) {
            
            double dt = t-out.r[nr].t0;
			double rem = dt/out.r[nr].dt;
			rem -= round(rem);double *data;
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
                case TXX:
                    data = tau_xx;
                    break;
                case TZZ:
                    data = tau_zz;
                    break;
                case TXZ:
                    data = tau_xz;
                    break;
                default:
                    break;
            }
            if ( ( fabs(dt) < 0.99*c.dt ) ||
                ( ( t > out.r[nr].t0 && out.r[nr].dt > 0.0 && ( fabs(rem) < 1.e-10 ) ) ) ) {
                if ( out.r[nr].type == TRACE )
                    write_trace(data, t, &out, nr);
                else if ( out.r[nr].type == SNAPSHOT ) {
//                    write_snapshot(data, t, &g, &out, nr, params.plotStrips);
                    write_snapshot_nc(data, t, &g, &out, nr, params.plotStrips);
                }
            }
        }
    }
	
	for ( size_t n=0; n<out.nrec; ++n )
        if ( out.r[n].type == TRACE )
            fclose( out.r[n].fid );
	
	for (size_t n=0; n<src.nsrc; ++n) free( src.s[n].fct );
	free ( src.s );
    free_fftw_data(&fd);
	
    if ( params.iwipe ) {
        free ( g.bc.gobx );
        free ( g.bc.gobz );
    }
    free ( delta );
    free ( Wtmp );
    free ( Ws );
    free ( W );
    free ( coeff_q_3 );
    free ( coeff_q_1 );
    free ( coeff_v_3 );
    free ( coeff_v_1 );
    free ( c.m3 );
    free ( c.m1 );
    free ( c.nk3 );
    free ( c.nk1 );
    free ( c.rho_f_j );
    free ( c.rho_f_i );
    free ( c.rho_j );
    free ( c.rho_i );
    free ( c.tau_s );
    free ( c.varphi );
    free ( c.alpha3 );
    free ( c.alpha1 );
    free ( c.M );
    free ( c.epsilon );
	free ( c.c55 );
	free ( c.c33 );
	free ( c.c13 );
	free ( c.c11 );
	
    if ( verbose ) fprintf(stdout, "\nEnd of computation\n");
    
	return 0;
}