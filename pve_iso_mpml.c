/*
 *  pve_iso_mpml.c
 *
 *  Created by Bernard Giroux on 11-13-27.
 *
 *  PoroViscoElastic wave propagation in 2D isotropic media, for one attenuation
 *  mechanism, on a staggered grid
 *
 *  Bernard Giroux
 *  INRS-ETE
 *
 *
 *  Code specs:  - language: ANSI C99
 *               - compiled with intel compiler 11.1 on a mac running OSX 10.6
 *               - external libraries: fftw ( http://www.fftw.org ), 
 *                       netCDF ( http://www.unidata.ucar.edu/software/netcdf/ )
 *  
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
 
 @ARTICLE{martin10,
 author = {R. Martin and D. Komatitsch and S. D. Gedney and E. Bruthiaux},
 title = {A High-Order Time and Space Formulation of the Unsplit Perfectly
 Matched Layer for the Seismic Wave Equation Using Auxiliary Differential
 Equations (ADE-PML)},
 journal = {CMES: Computer Modeling in Engineering \& Sciences},
 year = {2010},
 volume = {56},
 pages = {17--42},
 doi = {10.3970/cmes.2010.056.017},
 owner = {giroux},
 timestamp = {2011.03.08}
 }
 note: the above paper has a sign error in eq 32: should read + before c_{1,i}
 and another error in eq 33: \theta at the numerator should be (1-\theta).
 *
 *
 */

/*
 
 ------------------------------------------------------------------------------
 
 Usage
 
 pve_iso [options] -p parameter_file.dat
 
 options are:
 -h Print this message
 -v verbose mode (type twice for increased verbosity)
 -c checkpoint_file.dat   Restart run using data in checkpoint_file.dat (-p argument ignored)
 
 Format of parameter file:
 
 One line per parameter, first word is the value of the parameter, next is a
 keyword comprised between a # and a comma, and a comment can be added after the
 comma, i.e we have on each line
 
 value     # keyword, optional comment 
 
 Example showing available keywords:
 --- beginning of file ---
 model.dat      # model, model file
 source.dat     # source, source file
 output.dat     # output, description of output records
 run1           # basename, common for all output files
 1.2            # time, length of time window in s (default is 1.2)
 0.1            # dt, time step in ms (default is 0.1)
 1              # abs, apply absorbing boundary (0 or 1, 0 by default)
 20             # nl, number of absorbing layers (20 by default)
 1              # plstr, plot absorbing strips in snapshots (0 or 1, 0 by default)
 7              # kappa_max, PML parameter (1 by default)
 1              # alpha_max, PML parameter.  If 0, alpha = 1 everywhere in the PML; if 1 alpha_max = pi*f0
 0.001          # Rc, theoretical coefficient of reflection of PML layer (1e-5 by default)
 1              # saveEnergy, save E as fct of time
 0/500/0/500    # EnergyROI, region of interest for saving E
 1000           # chkpt_inc, checkpoint saving increment (0 by default, no checkpoint saving)
 --- end of file ---
 
 
 
 ------------------------------------------------------------------------------
 
 Format of model file
 
 --- beginning of file ---
 2 2            <-  nx and nz, number of nodes in x & z
 1.0 1.0        <-  dx and dz, grid step size in x & z, in meters
 0.0 0.0        <-  origin of the grid, in meters
 isotropic      <-  keyword for type of model (isotropic, anisotropic_vti)
 K_m(1,1) K_s(1,1) K_f(1,1) phi(1,1) mu(1,1) rho_s(1,1) rho_f(1,1) T(1,1) eta(1,1) kappa(1,1) Q(1,1) f0(1,1)
 K_m(1,2) K_s(1,2) K_f(1,2) phi(1,2) mu(1,2) rho_s(1,2) rho_f(1,2) T(1,2) eta(1,2) kappa(1,2) Q(1,2) f0(1,2)
 K_m(2,1) K_s(2,1) K_f(2,1) phi(2,1) mu(2,1) rho_s(2,1) rho_f(2,1) T(2,1) eta(2,1) kappa(2,1) Q(2,1) f0(2,1)
 K_m(2,2) K_s(2,2) K_f(2,2) phi(2,2) mu(2,2) rho_s(2,2) rho_f(2,2) T(2,2) eta(2,2) kappa(2,2) Q(2,2) f0(2,2)
 ...
 K_m(nx,nz) K_s(nx,nz) K_f(nx,nz) phi(nx,nz) mu(nx,nz) rho_s(nx,nz) rho_f(nx,nz) T(nx,nz) eta(nx,nz) kappa(nx,nz) Q(nx,nz) f0(nx,nz)
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
 f0        relaxation frequency                       [ Hz ]
 
 note: 1 cP = 1e-3 Pa.s ; 1 mD = 9.869233e-10 m^2
 
 *** Worthy option ***
 
 Layered models can be input by giving values for just a vertical profile; they 
 will be duplicated automatically.  For example:
 
 --- beginning of file ---
 2 2            <-  nx and nz, number of nodes in x & z
 1.0 1.0        <-  dx and dz, grid step size in x & z, in meters
 0.0 0.0        <-  origin of the grid, in meters
 isotropic      <-  keyword for type of model (isotropic, anisotropic_vti)
 K_m(1,1) K_s(1,1) K_f(1,1) phi(1,1) mu(1,1) rho_s(1,1) rho_f(1,1) T(1,1) eta(1,1) kappa(1,1) Q(1,1) f0(1,1)
 K_m(1,2) K_s(1,2) K_f(1,2) phi(1,2) mu(1,2) rho_s(1,2) rho_f(1,2) T(1,2) eta(1,2) kappa(1,2) Q(1,2) f0(1,2)
 ...
 K_m(1,nz) K_s(1,nz) K_f(1,nz) phi(1,nz) mu(1,nz) rho_s(1,nz) rho_f(1,nz) T(1,nz) eta(1,nz) kappa(1,nz) Q(1,nz) f0(1,nz)
 --- end of file ---
 
 
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

#include <sys/stat.h>
#include <errno.h>
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io_utils.h"
#include "pml.h"
#include "propagate.h"
#include "structs.h"
#include "src.h"

int verbose=0;

int main (int argc, char *argv[]) {
    
    struct inputParams params;
    struct computationVariables c;
    struct fac_pml fp1, fp23, fp4;
	struct mem_pml mem[3];
    struct materialProperties mp;
    struct fftw_data fd;
	struct sourceParams src;
    struct outputParams out;
	struct saveEnergy se;
	
    // -------------------------------------------------------------------------
    //
    // variables
    //
    // -------------------------------------------------------------------------
    
    const double pi = 4.0*atan(1.0);
//    const double L = 1.0;  // number of attenuation mechanisms
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
    
//    double *e;         // memory variable (one att. mechanism),
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
    
	size_t itstart = 0;
    
	params.ab = &g.ab;
    set_defaults(&g, &params);
    process_args(argc, argv, &params);
	
	if ( params.checkpoint == 1 ) {
		if ( verbose ) {
			fprintf(stdout, "\n *** pve_vti - PoroViscoElastic wave propagation in 2D isotropic media ***\n\n");
			fprintf(stdout, "  Starting with checkpoint file %s\n\n", params.checkpointfile);
		}
		read_checkpoint1(params.checkpointfile, &itstart, &g, &params);
		itstart++;
		if ( verbose ) {
			if ( verbose > 1 ) {
				fprintf(stdout, "  Model file: \t%s\n", params.modelfile);
				fprintf(stdout, "  Source file:\t%s\n", params.sourcefile);
				fprintf(stdout, "  Output file:\t%s\n", params.outputfile);
				fprintf(stdout, "  Time window:\t%lg (s)\n", params.duration);
				fprintf(stdout, "  Time step:  \t%lg (ms)\n", params.dt * 1.e3);
				fprintf(stdout, "  Model parameters: \t\tX\t\tZ\n");
				fprintf(stdout, "        grid points \t\t%zd\t\t%zd\n", g.nx, g.nz);
				fprintf(stdout, "        grid spacing\t\t%lg\t\t%lg\n", g.dx, g.dz);
				fprintf(stdout, "        grid origin \t\t%lg\t\t%lg\n", g.x0, g.z0);
				fprintf(stdout, "        pml padding\t\t%zd\t\t%zd\n", g.ab.np, g.ab.np);
				fprintf(stdout, "        pml - Rc       \t\t%lg\n", g.ab.Rc);
				fprintf(stdout, "        pml - kappa_max\t\t%lg\n", g.ab.kappa_max);
			}
			fprintf(stdout, "Allocating memory ... ");
			fflush(stdout);
		}			
	} else {
		
		strcpy(out.basename, params.basename);
		
		if ( verbose ) {
			fprintf(stdout, "\n *** pve_iso - PoroViscoElastic wave propagation in 2D isotropic media ***\n\n");
			if ( verbose > 1 ) {
				fprintf(stdout, "  Model file: \t%s\n", params.modelfile);
				fprintf(stdout, "  Source file:\t%s\n", params.sourcefile);
				fprintf(stdout, "  Output file:\t%s\n", params.outputfile);
				fprintf(stdout, "  Time window:\t%lg (s)\n", params.duration);
				fprintf(stdout, "  Time step:  \t%lg (ms)\n", params.dt * 1.e3);
				if ( params.chkpt_inc > 0 )
					fprintf(stdout, "  Checkpointing every %d iterations applied\n", params.chkpt_inc );
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
				fprintf(stdout, "        pml padding\t\t%zd\t\t%zd\n", g.ab.np, g.ab.np);
				fprintf(stdout, "        pml - Rc       \t\t%lg\n", g.ab.Rc);
				fprintf(stdout, "        pml - kappa_max\t\t%lg\n", g.ab.kappa_max);
				if ( g.ab.median ) fprintf(stdout, "        pml - median smoothing applied\n");
			}
			fprintf(stdout, "Allocating memory ... ");
			fflush(stdout);
		}
		
	}
	
	c.dt = params.dt;
	
    //
    //  Memory allocation
    //
    
    size_t nnodes = g.nx2 * g.nz2;
    
    if ( NULL == ( mp.K_m    = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.K_s    = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.K_f    = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.phi    = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.mu     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.rho_s  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.rho_f  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.T      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.eta    = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.kappa  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.Q      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.f0     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
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
	
	if ( NULL == ( W      = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( Ws     = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( Wtmp   = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( delta  = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	
	alloc_cpml(&fp1, &fp23, &fp4, mem, &g);
    c.mu = &(mp.mu[0]);
    
	if ( params.checkpoint == 1 ) {
		
		if ( verbose ) {
			fprintf(stdout, "done.\nReading main variables from checkpoint file ... ");
			fflush(stdout);
		}
		read_checkpoint2(params.checkpointfile, &g, &params, &src, &out, &se,
						 &fp1, &fp23, &fp4, mem, &c, coeff_v_i, coeff_v_j,
						 coeff_q_i, coeff_q_j, W, Ws, Wtmp, delta);
		
		if ( verbose ) {
			fprintf(stdout, "done.\n");
			fflush(stdout);
		}
	} else {
		if ( params.chkpt_inc > 0 ) params.checkpoint = 1;
		
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
		compute_src_fct(&src, params.dt, 0.5);
		normalize_src_fct(&src, &g);
		
		//
		// Read in output params
		//
		if ( verbose ) {
			fprintf(stdout, "done.\n");
			if ( verbose > 1 ) {
				fprintf(stdout, "Source has %zd component(s)\n", src.nsrc);
				double fac =1.;
				if ( src.nTemplate == 9 ) {
					fac = 2.;
					fprintf(stdout, "  Each component uses a %zd points template\n", src.nTemplate);
				}
				for ( size_t ns=0; ns<src.nsrc*src.nTemplate; ns+=src.nTemplate ) {
                    size_t nos = 1+(ns/src.nTemplate);
					fprintf(stdout, "    %zd - type:       \t%s\n", nos, typeSrc[src.s[ns].type]);
					fprintf(stdout, "    %zd - frequency:  \t%lg Hz\n", nos, src.s[ns].f);
					fprintf(stdout, "    %zd - strength:   \t%lg MPa\n", nos, src.s[ns].A*fac);
					fprintf(stdout, "    %zd - coordinates:\t(%lg, %lg)\n", nos, src.s[ns].x, src.s[ns].z);
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
						fprintf(stdout, "    %zd - Trace of %s at (%lg, %lg), ",
								nr+1, component[out.r[nr].comp], out.r[nr].x, out.r[nr].z);
						fprintf(stdout, "starting at %lg s & sampled at %lg ms\n",
								out.r[nr].t0, 1.e3*out.r[nr].dt);
					}
					else if ( out.r[nr].type == SNAPSHOT ) {
						if ( out.r[nr].dt <= 0.0 )
							fprintf(stdout, "    %zd - Snapshot of %s at time %lg s\n",
									nr+1, component[out.r[nr].comp], out.r[nr].t0);
						else
							fprintf(stdout, "    %zd - Snapshot of %s, starting at %lg s & sampled at %lg ms\n",
									nr+1, component[out.r[nr].comp], out.r[nr].t0, 1.e3*out.r[nr].dt);
					}
				}
			}
		}
		if ( params.saveEnergy == 1 ) {
			struct stat sb;
			if ( stat("E", &sb) < 0 ) {
				if ( errno == ENOENT ) {
					mkdir("E", S_IRWXU);
				}
			}
			char fname[80];
			sprintf(fname, "E/%s_E.dat", params.basename );
			se.fid = fopen( fname, "w+" );
			if ( verbose ) {
				fprintf(stdout, "Saving kinetic energy as fct of time in %s\n", fname);
				if (params.roiExmin<g.x0 || params.roiExmax>(g.x0+(g.nx-1)*g.dx) ||
					params.roiEzmin<g.z0 || params.roiEzmax>(g.z0+(g.nz-1)*g.dz)) {
					fprintf(stderr, "Region of interest for saving energy not fitting in the modeling grid\n");
					abort();
				}
				
				fprintf(stdout, "  Region of interest is %lg/%lg/%lg/%lg\n",
						params.roiExmin, params.roiExmax, params.roiEzmin, params.roiEzmax);
			}
			se.i1E = lround( (params.roiExmin-g.x0)/g.dx );
			se.i2E = lround( (params.roiExmax-g.x0)/g.dx );
			se.j1E = lround( (params.roiEzmin-g.z0)/g.dz );
			se.j2E = lround( (params.roiEzmax-g.z0)/g.dz );
		}
		
		if ( verbose ) fprintf(stdout, "Computing CPML parameters ... ");
		double alpha_max = pi*src.s[0].f;
		if ( g.ab.alpha_max == 0 ) alpha_max = 0.0;
		compute_cpml(&fp1, &fp23, &fp4, &g, params.dt, alpha_max, params.iwipe);
		
		if ( verbose ) {
			fprintf(stdout, "done.\nPreprocessing data ... ");
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
		// Compute elastic & other physical params
		//
		for ( size_t n=0; n<nnodes; ++n ) {
			double tau0 = 1./(2. * pi * mp.f0[n]);
			c.E[n] = mp.K_m[n] + 4./3.*mp.mu[n];                                 // eq (6)
			double D = mp.K_s[n] * (1. + mp.phi[n]*(mp.K_s[n]/mp.K_f[n] - 1. )); // eq (8)
			c.M[n] = mp.K_s[n]*mp.K_s[n] / (D - mp.K_m[n]);                      // eq (7)
			c.alpha[n] = 1. - mp.K_m[n]/mp.K_s[n];                               // eq (9)
			
			m[n] = mp.T[n]*mp.rho_f[n] / mp.phi[n];
			
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
		
		// smoothing absorbing region with median filter of increasing size (3x3 up to 9x9)
		if ( g.ab.median == 1 ) {
			
			pml_median_average( c.mu, &g );
			pml_median_average( c.E, &g );
			pml_median_average( c.M, &g );
			pml_median_average( c.alpha, &g );
			pml_median_average( c.varphi, &g );
			pml_median_average( c.tau_s, &g );
			pml_median_average( rho, &g );
			pml_median_average( mp.rho_f, &g );
			pml_median_average( mp.eta, &g );
			pml_median_average( mp.kappa, &g );
			pml_median_average( m, &g );
		}
		
		// check phase velocities and attenuation
		double c_p_m[6] = { 0., 0., 0., 1.e9, 1.e9, 1.e9 };
		double c_p_p[6] = { 0., 0., 0., 1.e9, 1.e9, 1.e9 };
		double c_s[6]   = { 0., 0., 0., 1.e9, 1.e9, 1.e9 };
		double a_p_m[6] = { 0., 0., 0., 1.e9, 1.e9, 1.e9 };
		double a_p_p[6] = { 0., 0., 0., 1.e9, 1.e9, 1.e9 };
		double a_s[6]   = { 0., 0., 0., 1.e9, 1.e9, 1.e9 };
		if ( verbose > 1 ) {
			for ( size_t n=0; n<nnodes; ++n ) {
				double tau0 = 1./(2. * pi * mp.f0[n]);
				double tau_e = tau0;
				if ( mp.Q[n] < 1.e4 ) {
					tau_e = tau0/mp.Q[n] * (sqrt(mp.Q[n]*mp.Q[n] + 1.) + 1.);
				}
				
				double f = 1.0;   // low frequency, 1 Hz
				double omega = 2.*pi*f;
				
				complex double rho_bar = mp.T[n]/mp.phi[n]*mp.rho_f[n] - I/omega * (mp.eta[n]/mp.kappa[n]);
				complex double rho_c = rho[n] - mp.rho_f[n]*mp.rho_f[n]/rho_bar;
				complex double M_c = c.M[n] / ( 1. + c.varphi[n] ) * (1. + I*omega*tau_e)/(1. + I*omega*c.tau_s[n]);
				complex double A = M_c*(rho[n]-2.*c.alpha[n]*mp.rho_f[n]) +
				rho_bar*(c.E[n] + c.alpha[n]*c.alpha[n]*M_c);
				
				complex double tmp = csqrt( A*A - 4.*M_c*c.E[n]*rho_c*rho_bar);
				complex double V_p_m = csqrt( ( A - tmp ) / (2.*rho_c*rho_bar) );
				complex double V_p_p = csqrt( ( A + tmp ) / (2.*rho_c*rho_bar) );
				complex double V_s = csqrt( mp.mu[n] / rho_c );
				
				double t = 1. / ( creal(1./V_p_m) );
				c_p_m[0] = t>c_p_m[0] ? t : c_p_m[0];
				c_p_m[3] = t<c_p_m[3] ? t : c_p_m[3];
				t = 1. / ( creal(1./V_p_p) );
				c_p_p[0] = t>c_p_p[0] ? t : c_p_p[0];
				c_p_p[3] = t<c_p_p[3] ? t : c_p_p[3];
				t = 1. / ( creal(1./V_s) );
				c_s[0]   = t>c_s[0] ? t : c_s[0];
				c_s[3]   = t<c_s[3] ? t : c_s[3];
				
				t = 17.372*pi*cimag(V_p_m)/creal(V_p_m);
				a_p_m[0] = t>a_p_m[0] ? t : a_p_m[0];
				a_p_m[3] = t<a_p_m[3] ? t : a_p_m[3];
				t = 17.372*pi*cimag(V_p_p)/creal(V_p_p);
				a_p_p[0] = t>a_p_p[0] ? t : a_p_p[0];
				a_p_p[3] = t<a_p_p[3] ? t : a_p_p[3];
				t = 17.372*pi*cimag(V_s)/creal(V_s);
				a_s[0]   = t>a_s[0] ? t : a_s[0];
				a_s[3]   = t<a_s[3] ? t : a_s[3];
				
				f = f0;   // src frequency
				omega = 2.*pi*f;
				
				rho_bar = mp.T[n]/mp.phi[n]*mp.rho_f[n] - I/omega * mp.eta[n]/mp.kappa[n];
				rho_c = rho[n] - mp.rho_f[n]*mp.rho_f[n]/rho_bar;
				M_c = c.M[n] / ( 1. + c.varphi[n] ) * (1. + I*omega*tau_e)/(1. + I*omega*c.tau_s[n]);
				A = M_c*(rho[n]-2.*c.alpha[n]*mp.rho_f[n]) +
				rho_bar*(c.E[n] + c.alpha[n]*c.alpha[n]*M_c);
				
				tmp = csqrt( A*A - 4.*M_c*c.E[n]*rho_c*rho_bar);
				V_p_m = csqrt( ( A - tmp ) / (2.*rho_c*rho_bar) );
				V_p_p = csqrt( ( A + tmp ) / (2.*rho_c*rho_bar) );
				V_s = csqrt( mp.mu[n] / rho_c );
				
				t = 1. / ( creal(1./V_p_m) );
				c_p_m[1] = t>c_p_m[1] ? t : c_p_m[1];
				c_p_m[4] = t<c_p_m[4] ? t : c_p_m[4];
				t = 1. / ( creal(1./V_p_p) );
				c_p_p[1] = t>c_p_p[1] ? t : c_p_p[1];
				c_p_p[4] = t<c_p_p[4] ? t : c_p_p[4];
				t = 1. / ( creal(1./V_s) );
				c_s[1]   = t>c_s[1] ? t : c_s[1];
				c_s[4]   = t<c_s[4] ? t : c_s[4];
				
				t = 17.372*pi*cimag(V_p_m)/creal(V_p_m);
				a_p_m[1] = t>a_p_m[1] ? t : a_p_m[1];
				a_p_m[4] = t<a_p_m[4] ? t : a_p_m[4];
				t = 17.372*pi*cimag(V_p_p)/creal(V_p_p);
				a_p_p[1] = t>a_p_p[1] ? t : a_p_p[1];
				a_p_p[4] = t<a_p_p[4] ? t : a_p_p[4];
				t = 17.372*pi*cimag(V_s)/creal(V_s);
				a_s[1]   = t>a_s[1] ? t : a_s[1];
				a_s[4]   = t<a_s[4] ? t : a_s[4];
				
				f = 1./c.dt;   // high frequency
				omega = 2.*pi*f;
				
				rho_bar = mp.T[n]/mp.phi[n]*mp.rho_f[n] - I/omega * mp.eta[n]/mp.kappa[n];
				rho_c = rho[n] - mp.rho_f[n]*mp.rho_f[n]/rho_bar;
				M_c = c.M[n] / ( 1. + c.varphi[n] ) * (1. + I*omega*tau_e)/(1. + I*omega*c.tau_s[n]);
				A = M_c*(rho[n]-2.*c.alpha[n]*mp.rho_f[n]) +
				rho_bar*(c.E[n] + c.alpha[n]*c.alpha[n]*M_c);
				
				tmp = csqrt( A*A - 4.*M_c*c.E[n]*rho_c*rho_bar);
				V_p_m = csqrt( ( A - tmp ) / (2.*rho_c*rho_bar) );
				V_p_p = csqrt( ( A + tmp ) / (2.*rho_c*rho_bar) );
				V_s = csqrt( mp.mu[n] / rho_c );
				
				t = 1. / ( creal(1./V_p_m) );
				c_p_m[2] = t>c_p_m[2] ? t : c_p_m[2];
				c_p_m[5] = t<c_p_m[5] ? t : c_p_m[5];
				t = 1. / ( creal(1./V_p_p) );
				c_p_p[2] = t>c_p_p[2] ? t : c_p_p[2];
				c_p_p[5] = t<c_p_p[5] ? t : c_p_p[5];
				t = 1. / ( creal(1./V_s) );
				c_s[2]   = t>c_s[2] ? t : c_s[2];
				c_s[5]   = t<c_s[5] ? t : c_s[5];
				
				t = 17.372*pi*cimag(V_p_m)/creal(V_p_m);
				a_p_m[2] = t>a_p_m[2] ? t : a_p_m[2];
				a_p_m[5] = t<a_p_m[5] ? t : a_p_m[5];
				t = 17.372*pi*cimag(V_p_p)/creal(V_p_p);
				a_p_p[2] = t>a_p_p[2] ? t : a_p_p[2];
				a_p_p[5] = t<a_p_p[5] ? t : a_p_p[5];
				t = 17.372*pi*cimag(V_s)/creal(V_s);
				a_s[2]   = t>a_s[2] ? t : a_s[2];
				a_s[5]   = t<a_s[5] ? t : a_s[5];
			}
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
		
		if ( params.saveEnergy == 1 ) {
			size_t npts = (se.i2E-se.i1E+1)*(se.j2E-se.j1E+1);
			if ( NULL == ( se.rho11 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
			if ( NULL == ( se.rho12 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
			if ( NULL == ( se.rho22 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
			// we should use averaged values here...
			for ( size_t i=se.i1E, n=0; i<=se.i2E; ++i ) {
				for ( size_t j=se.j1E; j<=se.j2E; ++j, ++n ) {
					size_t ij = (i+g.ab.np)*g.nz2+j+g.ab.np;
					se.rho12[n] = -mp.phi[ij]*mp.rho_f[ij]*1.e6*(mp.T[ij]-1.);  // eq. 7.160, Carcione (2007)
					se.rho11[n] = (1.-mp.phi[ij])*mp.rho_s[ij]*1.e6 - se.rho12[n]; // eq. 7.153
					se.rho22[n] = mp.phi[ij]*mp.rho_f[ij]*1.e6 - se.rho12[n];      // eq. 7.154
				}
			}
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
	    
		free ( mp.f0 );
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
			fprintf(stdout, "done.\n");
			if ( verbose > 1 ) {
				double f = 1./c.dt;   // high frequency
				fprintf(stdout, "Max/min phase velocities:\n");
				fprintf(stdout, "   Low frequency (1 Hz)    - slow P : %lf    %lf\n", c_p_m[0], c_p_m[3]);
				fprintf(stdout, "                             fast P : %lf    %lf\n", c_p_p[0], c_p_p[3]);
				fprintf(stdout, "                                  S : %lf    %lf\n", c_s[0], c_s[3]);
				fprintf(stdout, "   Src frequency (%lg Hz)   - slow P : %lf    %lf\n", f0, c_p_m[1], c_p_m[4]);
				fprintf(stdout, "                             fast P : %lf    %lf\n", c_p_p[1], c_p_p[4]);
				fprintf(stdout, "                                  S : %lf    %lf\n", c_s[1], c_s[4]);
				fprintf(stdout, "  High frequency (%lg Hz) - slow P : %lf    %lf\n", f, c_p_m[2], c_p_m[5]);
				fprintf(stdout, "                             fast P : %lf    %lf\n", c_p_p[2], c_p_p[5]);
				fprintf(stdout, "                                  S : %lf    %lf\n", c_s[2], c_s[5]);
				fprintf(stdout, "Max/min attenuation coefficient:\n");
				fprintf(stdout, "   Low frequency (1 Hz)    - slow P : %lf    %lf\n", a_p_m[0], a_p_m[3]);
				fprintf(stdout, "                             fast P : %lf    %lf\n", a_p_p[0], a_p_p[3]);
				fprintf(stdout, "                                  S : %lf    %lf\n", a_s[0], a_s[3]);
				fprintf(stdout, "   Src frequency (%lg Hz)   - slow P : %lf    %lf\n", f0, a_p_m[1], a_p_m[4]);
				fprintf(stdout, "                             fast P : %lf    %lf\n", a_p_p[1], a_p_p[4]);
				fprintf(stdout, "                                  S : %lf    %lf\n", a_s[1], a_s[4]);
				fprintf(stdout, "  High frequency (%lg Hz) - slow P : %lf    %lf\n", f, a_p_m[2], a_p_m[5]);
				fprintf(stdout, "                             fast P : %lf    %lf\n", a_p_p[2], a_p_p[5]);
				fprintf(stdout, "                                  S : %lf    %lf\n", a_s[2], a_s[5]);
			}
			fflush(stdout);
		}
		
		for ( size_t n=0; n<9*nnodes; ++n )
			W[n] = Ws[n] = delta[n] = 0.0;
	}
	
	
    tau_xx = &(W[0]);
    tau_zz = &(W[nnodes]);
    tau_xz = &(W[2*nnodes]);
    p      = &(W[3*nnodes]);
    v_x    = &(W[4*nnodes]);
    v_z    = &(W[5*nnodes]);
    q_x    = &(W[6*nnodes]);
    q_z    = &(W[7*nnodes]);
    
    vs_x    = &(Ws[4*nnodes]);
    vs_z    = &(Ws[5*nnodes]);
    qs_x    = &(Ws[6*nnodes]);
    qs_z    = &(Ws[7*nnodes]);
	
    tau_xxtmp = &(Wtmp[0]);
    tau_zztmp = &(Wtmp[nnodes]);
    tau_xztmp = &(Wtmp[2*nnodes]);
    ptmp      = &(Wtmp[3*nnodes]);
	
	/*
	 write_snapshot_nc(c.E, 0.0, &g, &out, 1000, params.plotStrips);
	 write_snapshot_nc(c.M, 0.0, &g, &out, 1001, params.plotStrips);
	 write_snapshot_nc(c.alpha, 0.0, &g, &out, 1002, params.plotStrips);
	 write_snapshot_nc(c.varphi, 0.0, &g, &out, 1003, params.plotStrips);
	 write_snapshot_nc(c.tau_s, 0.0, &g, &out, 1004, params.plotStrips);
	 */
	//write_snapshot_nc(c.rho_i, 0.0, &g, &out, 1005, params.plotStrips);
	/*
	 write_snapshot_nc(c.rho_j, 0.0, &g, &out, 1006, params.plotStrips);
	 write_snapshot_nc(c.rho_f_i, 0.0, &g, &out, 1007, params.plotStrips);
	 write_snapshot_nc(c.rho_f_j, 0.0, &g, &out, 1008, params.plotStrips);
	 write_snapshot_nc(c.nk_i, 0.0, &g, &out, 1009, params.plotStrips);
	 write_snapshot_nc(c.nk_j, 0.0, &g, &out, 1010, params.plotStrips);
	 write_snapshot_nc(c.mu_ij, 0.0, &g, &out, 1011, params.plotStrips);
	 write_snapshot_nc(c.mu, 0.0, &g, &out, 1012, params.plotStrips);
	 write_snapshot_nc(c.m_i, 0.0, &g, &out, 1013, params.plotStrips);
	 write_snapshot_nc(c.m_j, 0.0, &g, &out, 1014, params.plotStrips);
	 */
	
    init_fftw_data(Wtmp, &fd, &g);
	
	//
	//  Main loop
	//
	
    
    nsteps = params.duration / c.dt;
    if ( nsteps*c.dt < params.duration ) nsteps++;
    
    if ( verbose ) {
        fprintf(stdout, "done.\nStarting main loop (%zd iterations)\n", nsteps);
        fflush(stdout);
    }
    for ( size_t it=itstart; it<nsteps; ++it ) {
        
		double t = c.dt*(it+1);
		if ( (it+1)%50 == 0 && verbose ) { fprintf(stdout, "  Iteration %zd, t = %g s\n", it+1, t); fflush(stdout); }
		
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
		
		// update PML memory variables
		for ( size_t n=0; n<(2*g.ab.np+1)*g.nz2; ++n ) {
			mem[0].dx_txx[n] = mem[1].dx_txx[n];
			mem[0].dx_p[n]   = mem[1].dx_p[n];
			mem[0].dx_txz[n] = mem[1].dx_txz[n];
			mem[0].dx_vx[n]  = mem[1].dx_vx[n];
			mem[0].dx_qx[n]  = mem[1].dx_qx[n];
			mem[0].dx_vz[n]  = mem[1].dx_vz[n];
		}
		for ( size_t n=0; n<(2*g.ab.np+1)*g.nx2; ++n ) {
			mem[0].dz_txz[n] = mem[1].dz_txz[n];
			mem[0].dz_tzz[n] = mem[1].dz_tzz[n];
			mem[0].dz_p[n]   = mem[1].dz_p[n];
			mem[0].dz_vz[n]  = mem[1].dz_vz[n];
			mem[0].dz_qz[n]  = mem[1].dz_qz[n];
			mem[0].dz_vx[n]  = mem[1].dz_vx[n]; 
		}			
        
        // eq (80) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n )
            Wtmp[n] = Ws[n];  // -> using Wtmp rather that Ws for first ride because fftw plans use Wtmp
        
        propagateCPML(Wtmp, delta, &g, &c, &fd, mem, &fp1);
        
        // eq (81) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt/2. * delta[n];
			W[n]    = Ws[n] + c.dt/6. * delta[n];
		}
		// add source term
		add_src_rk4(&src, W, Wtmp, c.dt/6., c.dt/2., it, nnodes);
        
        propagateCPML(Wtmp, delta, &g, &c, &fd, mem, &fp23);
        
        // eq (82) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt/2. * delta[n];
			W[n]   +=         c.dt/3. * delta[n];
		}
		// add source term
		add_src_rk4(&src, W, Wtmp, c.dt/3., c.dt/2., it, nnodes);
        
        propagateCPML(Wtmp, delta, &g, &c, &fd, mem, &fp23);
        
        // eq (83) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt    * delta[n];
			W[n]   +=         c.dt/3. * delta[n];
		}
		// add source term
		add_src_rk4(&src, W, Wtmp, c.dt/3., c.dt, it, nnodes);
        
        propagateCPML(Wtmp, delta, &g, &c, &fd, mem, &fp4);
        
        // eq (79) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n )
            W[n] += c.dt/6. * delta[n];
		// add source term
		add_src_rk4(&src, W, Wtmp, c.dt/6., 0.0, it, nnodes);
		
		for ( size_t nr=0; nr<out.nrec; ++nr ) {
            
			double dt = t-out.r[nr].t0+1.e-15;
            double rem = fmod(1.e6*dt, 1.e6*out.r[nr].dt);
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
                case TXX:
                    data = tau_xx;
                    break;
                case TZZ:
                    data = tau_zz;
                    break;
                case TXZ:
                    data = tau_xz;
                    break;
                case P:
                    data = p;
                    break;
                default:
					data = NULL;
                    break;
            }
			if ( data == NULL ) continue;
			
            if ( ( fabs(dt) < 0.99*c.dt ) ||
                ( ( t > out.r[nr].t0 && out.r[nr].dt > 0.0 && ( fabs(rem) < 1.e-8 ) ) ) ) {
				
                if ( out.r[nr].type == TRACE )
                    write_trace(data, t, &out, nr);
                else if ( out.r[nr].type == SNAPSHOT ) {
					if ( verbose > 1 )
                        printf("Writing snapshot of %s, t = %g s\n", component[out.r[nr].comp], t);
                    // write_snapshot(data, t, &g, &out, nr, params.plotStrips);
                    write_snapshot_nc(data, t, &g, &out, nr, params.plotStrips);
				}
            }
        }
		
		if ( params.saveEnergy ) {
			double E = 0.0;
			for ( size_t i=se.i1E, n=0; i<=se.i2E; ++i ) {
				for ( size_t j=se.j1E; j<=se.j2E; ++j, ++n ) {
					size_t ij = (i+g.ab.np)*g.nz2+j+g.ab.np;
					
					// eq. 7.145 of Carcione (2007)  book
					E += 0.5*(se.rho11[n]*( v_x[ij]*v_x[ij] + v_z[ij]*v_z[ij] ));
					E += (se.rho12[n]*( v_x[ij]*q_x[ij] + v_z[ij]*q_z[ij] ));
					E += 0.5*(se.rho22[n]*( q_x[ij]*q_x[ij] + q_z[ij]*q_z[ij] ));
				}
			}
			fprintf(se.fid, "%le   %le\n", t, E);
		}
		
		if ( params.checkpoint==1 && it>0 && it%params.chkpt_inc==0 ) {
			if ( verbose > 1 )
				printf("Saving checkpoint data at iteration %zd\n", it);
			save_checkpoint(it, &g, &params, &src, &out, &se, &fp1, &fp23, &fp4,
							mem, &c, coeff_v_i, coeff_v_j,
							coeff_q_i, coeff_q_j, W, Ws, Wtmp, delta);
		}
    }
	
	for ( size_t n=0; n<out.nrec; ++n )
        if ( out.r[n].type == TRACE )
            fclose( out.r[n].fid );
	
	if ( params.saveEnergy == 1 ) {
		fclose( se.fid );
		free ( se.rho11 );
		free ( se.rho12 );
		free ( se.rho22 );
	}
	
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
    
	free_cpml(&fp1, &fp23, &fp4, mem);
	
    if ( verbose ) fprintf(stdout, "\nEnd of computation\n");
    
    return 0;
}
