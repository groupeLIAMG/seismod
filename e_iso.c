/*
 *  e_iso.c
 *
 * Elastic wave propagation in 2D isotropic media, on a staggered grid 
 *
 *  Code specs:  - language: ANSI C99
 *               - compiled with intel compiler 11.1 on a mac running OSX 10.6
 *               - external libraries:
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
 *   Reference paper
 *
 @ARTICLE{collino01,
 author = {Francis Collino and Chrysoula Tsogka},
 title = {Application of the perfectly matched absorbing layer model to the
 linear elastodynamic problem in anisotropic heterogeneous media},
 journal = {Geophysics},
 year = {2001},
 volume = {66},
 pages = {294-307},
 number = {1},
 doi = {10.1190/1.1444908},
 publisher = {SEG},
 url = {http://link.aip.org/link/?GPY/66/294/1}
 }
 *
 *  Created by Bernard Giroux on 11-02-11.
 *  Copyright 2011 INRS ETE. All rights reserved.
 *
 */


/*
 ------------------------------------------------------------------------------
 
 Usage
 
 e_iso [options] -p parameter_file.dat
 
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
 20                 # nl, number of PML layers (20 by default)
 --- end of file ---

 
 ------------------------------------------------------------------------------
 
 Format of model file
 
 --- beginning of file ---
 2 2            <-  nx and nz, number of nodes in x & z
 1.0 1.0        <-  dx and dz, grid step size in x & z, in meters
 0.0 0.0        <-  origin of the grid, in meters
 isotropic_elastic      <-  keyword for type of model 
 Vp(1,1) Vs(1,1) rho(1,1)
 Vp(1,2) Vs(1,2) rho(1,2)
 Vp(2,1) Vs(2,1) rho(2,1)
 Vp(2,2) Vs(2,2) rho(2,2)
 ...
 Vp(nx,nz) Vs(nx,nz) rho(nx,nz)
 --- end of file ---
 
 Variable  Description                                 Units
 
 Vp        P-wave velocity                           [ m/s ]
 Vs        S-wave velocity                           [ m/s ]
 rho       density                                [ kg/m^3 ]
 
 *** Worthy option ***
 
 Layered models can be input by giving values for just a vertical profile; they 
 will be duplicated automatically.  For example:
 
 --- beginning of file ---
 2 2            <-  nx and nz, number of nodes in x & z
 1.0 1.0        <-  dx and dz, grid step size in x & z, in meters
 0.0 0.0        <-  origin of the grid, in meters
 isotropic_elastic      <-  keyword for type of model
 Vp(1,1) Vs(1,1) rho(1,1)
 Vp(1,2) Vs(1,2) rho(1,2)
 ...
 Vp(1,nz) Vs(1,nz) rho(1,nz)
 --- end of file ---
 
 
 ------------------------------------------------------------------------------
 
 
 */

char *parFormat =
"Format of parameter file:\n\
\n\
One line per parameter, first word is the value of the parameter,\n\
next is a keyword comprise between a # and a comma, and a comment\n\
can be added after the comma, i.e we have on each line\n\n\
value     # keyword, optional comment\n\n\
Example showing available keywords:\n\
--- beginning of file ---\n\
model.dat      # model, model file\n\
source.dat     # source, source file\n\
output.dat     # output, description of output records\n\
run1           # basename, common for all output files\n\
1.2            # time, length of time window in s (default is 1.2)\n\
0.1            # dt, time step in ms (default is 0.1)\n\
0              # segy, save traces in segy format (0 or 1, 0 by default)\n\
1              # shotpoint, Shotpoint number for segy file (1 by default)\n\
0              # simulateCMP, for 1D model save SEGY as a CMP\n\
1              # abs, apply absorbing boundary (0 or 1, 0 by default)\n\
20             # nl, number of absorbing layers (20 by default)\n\
1              # plstr, plot absorbing strips in snapshots (0 or 1, 0 by default)\n\
--- end of file ---\n";

char *modelFormat =
"Format of model file\n\
\n\
--- beginning of file ---\n\
2 2                  <-  nx and nz, number of nodes in x & z\n\
1.0 1.0              <-  dx and dz, grid step size in x & z, in meters\n\
0.0 0.0              <-  origin of the grid, in meters\n\
isotropic_elastic      <-  keyword for type of model (isotropic, anisotropic_vti)\n\
V_p(1,1) V_s(1,1) rho(1,1)\n\
V_p(1,2) V_s(1,2) rho(1,2)\n\
V_p(2,1) V_s(2,1) rho(2,1)\n\
V_p(2,2) V_s(2,2) rho(2,2)\n\
...\n\
V_p(nx,nz) V_s(nx,nz) rho(nx,nz)\n\
--- end of file ---\n\
\n\
Variable  Description                                 Units\n\
\n\
V_p       P-wave velocity                           [ m/s ]\n\
V_s       S-wave velocity                           [ m/s ]\n\
rho       density                                [ kg/m^3 ]\n\
\n\
*** Worthy option ***\n\
\n\
Layered models can be input by giving values for just a vertical profile; they\n\
will be duplicated automatically.  For example:\n\
\n\
--- beginning of file ---\n\
2 2            <-  nx and nz, number of nodes in x & z\n\
1.0 1.0        <-  dx and dz, grid step size in x & z, in meters\n\
0.0 0.0        <-  origin of the grid, in meters\n\
isotropic_elastic      <-  keyword for type of model (isotropic, anisotropic_vti)\n\
V_p(1,1) V_s(1,1) rho(1,1)\n\
V_p(1,2) V_s(1,2) rho(1,2)\n\
...\n\
V_p(1,nz) V_s(1,nz) rho(1,nz)\n\
--- end of file ---\n";

char *srcFormat =
"Format of source file\n\
\n\
--- beginning of file ---\n\
2              <- number of source points\n\
1              <- Number of pts in source template, can be 1 or 9 (Lin et Thylén, 2009)\n\
\n\
Sz             <- Source component of 1st source point\n\
10.0           <- Strength in MPa at 1st source point\n\
50.0           <- Frequency in Hz at 1st source point\n\
50.0 45.0      <- x z coordinates at 1st source point in meters\n\
\n\
Sx             <- Source component of 2nd source point\n\
1.0            <- Strength in MPa at 2nd source point\n\
50.0           <- Frequency in Hz at 2nd source point\n\
50.0 45.0      <- x z coordinates at 2nd source point in meters\n\
--- end of file ---\n\
\n\
Components can be Sx, Sz, Sxz or Bulk\n\
\n\
@ARTICLE{lin09,\n\
author = {Zhili Lin and Lars Thylén},\n\
title = {An analytical derivation of the optimum source patterns for the\n\
pseudospectral time-domain method},\n\
journal = {Journal of Computational Physics},\n\
year = {2009},\n\
volume = {228},\n\
pages = {7375 - 7387},\n\
number = {19},\n\
doi = {10.1016/j.jcp.2009.06.033}\n\
}\n";


char *rcvFormat =
"Format of file defining the outputs\n\
\n\
--- beginning of file ---\n\
2              <- number of records\n\
\n\
Vz             <- particle velocity component of 1st record\n\
Trace          <- type of 1st record\n\
10.0 20.0      <- X Z coordinates at 1st record point in meters\n\
0.0            <- start time of 1st recording, in seconds\n\
2.0            <- time sampling of 1st record, in miliseconds\n\
\n\
Vx             <- particle velocity component of 2nd record\n\
Snapshot       <- type of 2nd record\n\
0.0            <- start time of snapshot(s) (2nd recording), in seconds\n\
2.0            <- time sampling of snapshot(s) 2nd record, in miliseconds.\n\
Give 0.0 or less for just one snapshot at time given above.\n\
--- end of file ---\n\
\n\
Trace files are ascii files and snapshots are binary grid files in netCDF format\n\
(COARDS-compliant).\n\
\n\
Files names are\n\
\n\
basename_NO.trc            for traces, where NO is the record number\n\
basename_NO_TIME.nc        for snapshots, where NO is record number and TIME is\n\
the record time\n\
\n\
Trace files are two columns, first for time and second for velocity component.\n\
Velocity components in traces are _not_ interpolated at exact X and Z\n\
coordinates, but rather taken at the closest grid point.\n\
Snapshot  can be viewed in matlab, with Paraview ( http://www.paraview.org/ )\n\
or by using GMT routines ( http://www.soest.hawaii.edu/gmt/ ) to create plots.\n\
\n\
Units of particle velocity are mm/s, time in s and distances in m.\n";



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io_utils.h"
#include "pml.h"
#include "propagate.h"
#include "structs.h"
#include "src.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int verbose=0;

int main (int argc, char *argv[]) {
	
	struct inputParams params;
	struct sourceParams src;
    struct outputParams out;

    // -------------------------------------------------------------------------
    //
    // variables
    //
    // -------------------------------------------------------------------------
	
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
	double *vx, *vz, *txx, *tzz, *txz, *la, *mu, *l2m, *b, *bm;
    
	//
	// Material properties (input)
	//
	double *vp, *vs, *rho, *mu_tmp, *la_tmp, *l2m_tmp;
	
	//
    // pml variables
	//
    double *d;  // damping factor at grid nodes
    double *dh; // damping factor at half grid nodes at beggining of grid
	double *dH; // damping factor at half grid nodes at end of grid
    double *vx_x, *vx_z;  // 
    double *vz_x, *vz_z;
    double *txx_x, *txx_z;
    double *tzz_x, *tzz_z;
    double *txz_x, *txz_z;
	
	params.ab = &g.ab;
	set_defaults(&g, &params);
    process_args(argc, argv, &params);
    double dt = params.dt;
    size_t Npml = params.ab->np;
	g.ab.np = 0;  // domain not extended here, pml handled independently
    strcpy(out.basename, params.basename);

	if ( verbose ) {
        fprintf(stdout, "\n *** e_iso - Elastic wave propagation in 2D isotropic media ***\n\n");
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
	
	if ( g.dx != g.dz ) {
		fprintf(stderr, "Error: grid cell size must be equal in X and Z, aborting.\n");
		abort();
	}
	g.nx2 = g.nx;
	g.nz2 = g.nz;
	double h = g.dx;
	
	if ( verbose ) {
        fprintf(stdout, "done.\n");
        if ( verbose > 1 ) {
            fprintf(stdout, "  Model parameters: \t\tX\t\tZ\n");
            fprintf(stdout, "        grid points \t\t%zd\t\t%zd\n", g.nx, g.nz);
            fprintf(stdout, "        grid spacing\t\t%lg\t\t%lg\n", g.dx, g.dz);
            fprintf(stdout, "        grid origin \t\t%lg\t\t%lg\n", g.x0, g.z0);
            fprintf(stdout, "        grid padding\t\t%zd\t\t%zd\n", g.ab.np, g.ab.np);
        }
        fprintf(stdout, "Allocating memory ... ");
        fflush(stdout);
    }
	
	vx =  (double *)malloc(g.nx*g.nz*sizeof(double)); // particle velocities
    vz =  (double *)malloc(g.nx*g.nz*sizeof(double));
    txx = (double *)malloc(g.nx*g.nz*sizeof(double)); // stresses
    tzz = (double *)malloc(g.nx*g.nz*sizeof(double));
    txz = (double *)malloc(g.nx*g.nz*sizeof(double));
    la =  (double *)malloc(g.nx*g.nz*sizeof(double)); // Lamé cte (lambda)
    mu =  (double *)malloc(g.nx*g.nz*sizeof(double)); // Lamé cte
    l2m = (double *)malloc(g.nx*g.nz*sizeof(double)); // lambda+2*mu
    b =   (double *)malloc(g.nx*g.nz*sizeof(double)); // 1/rho
    bm =  (double *)malloc(g.nx*g.nz*sizeof(double)); // 1/rho (moyenné)

	la_tmp = txx;
	mu_tmp = tzz;
	l2m_tmp = txz;
	vp = vx;
	vs = vz;
	rho = bm;
    
    vx_x  = (double *)malloc((4*Npml*Npml+2*Npml*g.nx+2*Npml*g.nz)*sizeof(double));
    vx_z  = (double *)malloc((4*Npml*Npml+2*Npml*g.nx+2*Npml*g.nz)*sizeof(double));
    vz_x  = (double *)malloc((4*Npml*Npml+2*Npml*g.nx+2*Npml*g.nz)*sizeof(double));
    vz_z  = (double *)malloc((4*Npml*Npml+2*Npml*g.nx+2*Npml*g.nz)*sizeof(double));
    txx_x = (double *)malloc((4*Npml*Npml+2*Npml*g.nx+2*Npml*g.nz)*sizeof(double));
    txx_z = (double *)malloc((4*Npml*Npml+2*Npml*g.nx+2*Npml*g.nz)*sizeof(double));
    tzz_x = (double *)malloc((4*Npml*Npml+2*Npml*g.nx+2*Npml*g.nz)*sizeof(double));
    tzz_z = (double *)malloc((4*Npml*Npml+2*Npml*g.nx+2*Npml*g.nz)*sizeof(double));
    txz_x = (double *)malloc((4*Npml*Npml+2*Npml*g.nx+2*Npml*g.nz)*sizeof(double));
    txz_z = (double *)malloc((4*Npml*Npml+2*Npml*g.nx+2*Npml*g.nz)*sizeof(double));
	
	//
    // Read in model
    // 
    
    if ( verbose ) {
        fprintf(stdout, "done.\nReading properties of materials ... ");
        fflush(stdout);
    }
    read_model_e(params.modelfile, &g, vp, vs, rho);
	
	double vpmax = 0.0;
    for (size_t i=0; i<g.nx; ++i) {
        for (size_t j=0; j<g.nz; ++j) {
			size_t ind = i*g.nz+j;
            vpmax = vp[ind]>vpmax ? vp[ind] : vpmax;
            la_tmp[ind] = rho[ind]*(vp[ind]*vp[ind]-2.*vs[ind]*vs[ind]);
            mu_tmp[ind] = rho[ind]*vs[ind]*vs[ind];
            l2m_tmp[ind] = la_tmp[ind] + 2.*mu_tmp[ind];
            b[ind] = 1./rho[ind];
        }
    }
    
	//
	// averaging material properties at staggered nodes
	//
    for (size_t i=0; i<g.nx; ++i) {
        for (size_t j=0; j<g.nz-1; ++j) {
			mu[i*g.nz+j] = ( mu_tmp[i*g.nz+j] + mu_tmp[i*g.nz+j+1] )/2.0;
		}
		mu[i*g.nz+g.nz-1] = mu_tmp[i*g.nz+g.nz-1];
	}
	for (size_t j=0; j<g.nz; ++j) {
		for (size_t i=0; i<g.nx-1; ++i) {
			la[i*g.nz+j] = ( la_tmp[i*g.nz+j] + la_tmp[(i+1)*g.nz+j] )/2.0;
			l2m[i*g.nz+j] = ( l2m_tmp[i*g.nz+j] + l2m_tmp[(i+1)*g.nz+j] )/2.0;
		}
		la[(g.nx-1)*g.nz+j] = la_tmp[(g.nx-1)*g.nz+j];
		l2m[(g.nx-1)*g.nz+j] = l2m_tmp[(g.nx-1)*g.nz+j];
	}
    for (size_t i=0; i<g.nx-1; ++i) {
        for (size_t j=0; j<g.nz-1; ++j) {
			bm[i*g.nz+j] = ( b[i*g.nz+j]   + b[(i+1)*g.nz+j] + 
							b[i*g.nz+j+1] + b[(i+1)*g.nz+j+1] )/4.0;
		}
		bm[i*g.nz+g.nz-1] = ( b[i*g.nz+g.nz-1] + b[(i+1)*g.nz+g.nz-1] )/2.0;
	}
	for (size_t j=0; j<g.nz-1; ++j) {
		bm[(g.nx-1)*g.nz+j] = ( b[(g.nx-1)*g.nz+j] + b[(g.nx-1)*g.nz+j+1] )/2.0;
	}
	bm[g.nx*g.nz-1] = b[g.nx*g.nz-1];
	
    
//    FILE *fid = fopen("model.out", "wt");
//    for (size_t i=0; i<g.nx; ++i) {
//        for (size_t j=0; j<g.nz; ++j) {
//            size_t ind = i*g.nz+j;
//            fprintf(fid, "%g  %g  %g  %g  %g  %g\n", g.x0+i*g.dx, g.z0+j*g.dz, la[ind], mu[ind], l2m[ind], bm[ind]);
//        }
//    }
//    fclose(fid);

    

	if ( dt > h / ( sqrt(2.0) * vpmax ) ) {
		fprintf(stderr, "\n *** Warning!  Time step too large, aborting ***\n");
		abort();
	}
    
    if ( verbose ) fprintf(stdout, "done.\nApplying PML absorbing boundaries (%zd layers)", Npml);
    d  = (double *)malloc(Npml*sizeof(double));
    dh = (double *)malloc(Npml*sizeof(double));
    dH = (double *)malloc((Npml+1)*sizeof(double));
		
    vpmax = 0.0;
    for ( size_t n=0; n<g.nx*g.nz; ++n ) vpmax += vp[n];
    vpmax /= (double)g.nx*g.nz;
		
    double R = 0.001;
    int pow_d = 2;
    double delta = Npml*h;
    double d0 = log(1./R) * 3./2. *vpmax/delta;
    for ( size_t n=0; n<Npml; ++n ) {
        d[n] = d0*pow( ((Npml-n)*h/delta), pow_d);
        dh[n] = d0*pow( (((Npml-n)*h-h/2.)/delta), pow_d);
    }
    for ( size_t n=0; n<Npml+1; ++n ) {
        dH[n] = d0*pow( (((Npml-n)*h+h/2.)/delta), pow_d);
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
	compute_src_fct(&src, params.dt, 1.0);
	
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
            for ( size_t ns=0; ns<src.nsrc; ++ns ) {
                fprintf(stdout, "    %zd - type:       \t%s\n", ns+1, typeSrc[src.s[ns].type]);
                fprintf(stdout, "    %zd - frequency:  \t%lg Hz\n", ns+1, src.s[ns].f);
                fprintf(stdout, "    %zd - strength:   \t%lg MPa\n", ns+1, src.s[ns].A*1.e3*fac);
                fprintf(stdout, "    %zd - coordinates:\t(%lg, %lg)\n", ns+1, src.s[ns].x, src.s[ns].z);
            }
        }
        fprintf(stdout, "Reading output parameters ... ");
        fflush(stdout);
    }
    read_output(params.outputfile, &g, &out);
    check_dt_trc(&out, dt);

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
        fflush(stdout);
    }
    
	for (size_t n=0; n<g.nx*g.nz; ++n) {
        vx[n] = vz[n] = txx[n] = tzz[n] = txz[n] = 0.0;
    }
	
	for (size_t n=0; n<4*Npml*Npml+(2*Npml*g.nx)+(2*Npml*g.nz); ++n) {
        vx_x[n] = vz_x[n] = txx_x[n] = tzz_x[n] = txz_x[n] = 0.0;
        vx_z[n] = vz_z[n] = txx_z[n] = tzz_z[n] = txz_z[n] = 0.0;
    }
    
	
    nsteps = params.duration / dt;
    if ( nsteps*dt < params.duration ) nsteps++;

#ifdef _OPENMP
    if ( verbose ) {
        fprintf(stdout, "Will be running threaded (OpenMP) version with 3 threads\n");
    }
    omp_set_num_threads( 3 );
#endif
    
    if ( verbose ) {
        fprintf(stdout, "Starting main loop (%zd iterations)\n", nsteps);
        fflush(stdout);
    }
	for (size_t it=0; it<nsteps; ++it) {
        
		double t = dt*(it+1);
		if ( (it+1)%50 == 0 && verbose ) { 
            fprintf(stdout, "  Iteration %zd, t = %g s\n", it+1, t); fflush(stdout);
        }
#pragma omp parallel default(shared)
        {
#pragma omp sections
            {
        // ------------------------------------------------------------
        // Vx
        // ------------------------------------------------------------
#pragma omp section
        update_vx(vx, txx, txz, b, txx_x, txx_z, txz_x, txz_z, &g, dt, h, Npml);
		
        // ------------------------------------------------------------
        // Vz
        // ------------------------------------------------------------
#pragma omp section
        update_vz(vz, txz, tzz, bm, txz_x, txz_z, tzz_x, tzz_z, &g, dt, h, Npml);
        
#pragma omp section
        pml_v(vx_x, vx_z, vz_x, vz_z, txx_x, txx_z, tzz_x, tzz_z, txz_x, txz_z,
              d, dh, dH, txx, tzz, txz, b, bm, g.nx, g.nz, Npml, dt, h);
        
            } /*-- End of sections block --*/
#pragma omp sections
            {
        // ------------------------------------------------------------
        // tau_xx and tau_zz
        // ------------------------------------------------------------
#pragma omp section
        update_txxzz(txx, tzz, vx, vz, la, l2m, vx_x, vx_z, vz_x, vz_z, &g, dt, h, Npml);
        
        // ------------------------------------------------------------
        // tau_xz
        // ------------------------------------------------------------
#pragma omp section
        update_txz(txz, vx, vz, mu, vx_x, vx_z, vz_x, vz_z, &g, dt, h, Npml);
        
#pragma omp section
        pml_t(vx_x, vx_z, vz_x, vz_z, txx_x, txx_z, tzz_x, tzz_z, txz_x, txz_z,
              d, dh, dH, vx, vz, la, l2m, mu, g.nx, g.nz, Npml, dt, h);
		
            } /*-- End of sections block --*/
        } /*-- End of parallel region --*/
        add_src(&src, txx, tzz, txz, dt, it);
			
        // ------------------------------------------------------------
        // write records
        // ------------------------------------------------------------
		for ( size_t nr=0; nr<out.nrec; ++nr ) {
            
			double dt2 = t-out.r[nr].t0+1.e-15;
            double rem = fmod(1.e6*dt2, 1.e6*out.r[nr].dt);
			double *data;
            switch (out.r[nr].comp) {
                case VX:
                    data = vx;
                    break;
                case VZ:
                    data = vz;
                    break;
                default:
                    break;
            }
			
            if ( ( fabs(dt2) < 0.99*dt ) ||
                ( ( t > out.r[nr].t0 && out.r[nr].dt > 0.0 && ( fabs(rem) < 1.e-8 ) ) ) ) {
				
                if ( out.r[nr].type == TRACE )
                    write_trace(data, t, &out, nr);
                else if ( out.r[nr].type == SNAPSHOT ) {
					if ( verbose > 1 )
                        printf("Writing snapshot of %s, t = %g s\n", component[out.r[nr].comp], t);
					//                    write_snapshot(data, t, &g, &out, nr, params.plotStrips);
                    write_snapshot_nc(data, t, &g, &out, nr, params.plotStrips);
				}
            }
		}			
	}	
	
	free( vx );
    free( vz );
    free( txx );
    free( tzz );
    free( txz );
    free( la );
    free( mu );
    free( l2m );
    free( b );
    free( bm );
	
    free( vx_x );
    free( vx_z );
    free( vz_x );
    free( vz_z );
    free( txx_x );
    free( txx_z );
    free( tzz_x );
    free( tzz_z );
    free( txz_x );
    free( txz_z );
    free( d );
    free( dh );
    free( dH );
	
    if ( verbose ) fprintf(stdout, "\nEnd of computation\n");    
    return 0;
}
