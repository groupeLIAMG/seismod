//
//  a_iso_3d - Acoustic wave propagation in 3D isotropic media, on a staggered grid
//
//  Created by Bernard Giroux on 2012-11-13.
//
//

/*
 * Copyright (c) 2012, Bernard Giroux
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
 ------------------------------------------------------------------------------
 
 Usage
 
 a_iso_3d [options] -p parameter_file.dat
 
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
 2 2 2              <-  nx,ny and nz, number of nodes in x, y & z
 1.0 1.0 1.0        <-  dx, dy and dz, grid step size in x, y & z, in meters
 0.0 0.0 0.0        <-  origin of the grid, in meters
 isotropic_acoustic_3d      <-  keyword for type of model
 Vp(1,1,1) rho(1,1,1)
 Vp(1,1,2) rho(1,1,2)
 Vp(1,2,1) rho(1,2,1)
 Vp(1,2,2) rho(1,2,2)
 ...
 Vp(nx,ny,nz) rho(nx,ny,nz)
 --- end of file ---
 
 Variable  Description                                 Units
 
 Vp        P-wave velocity                           [ m/s ]
 rho       density                                [ kg/m^3 ]
 
--- end of file ---
 
 
 ------------------------------------------------------------------------------
 
 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io_utils.h"
#include "structs.h"
#include "src.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int verbose=0;






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
1              # abs, apply absorbing boundary (0 or 1, 0 by default)\n\
20             # nl, number of PML layers (20 by default)\n\
1              # plstr, plot absorbing strips in snapshots (0 or 1, 0 by default)\n\
7              # kappa_max, PML parameter (1 by default)\n\
1              # alpha_max, PML parameter.  If 0, alpha = 1 everywhere in the PML; if 1 alpha_max = pi*f0\n\
0.001          # Rc, theoretical coeff of reflection of PML layer (1e-5 by default)\n\
--- end of file ---\n";

char *modelFormat =
" --- beginning of file ---\n\
2 2 2              <-  nx,ny and nz, number of nodes in x, y & z\n\
1.0 1.0 1.0        <-  dx, dy and dz, grid step size in x, y & z, in meters\n\
0.0 0.0 0.0        <-  origin of the grid, in meters\n\
isotropic_acoustic_3d      <-  keyword for type of model\n\
Vp(1,1,1) rho(1,1,1)\n\
Vp(1,1,2) rho(1,1,2)\n\
Vp(1,2,1) rho(1,2,1)\n\
Vp(1,2,2) rho(1,2,2)\n\
...\n\
Vp(nx,ny,nz) rho(nx,ny,nz)\n\
--- end of file ---\n\
\n\
Variable  Description                                 Units\n\
\n\
Vp        P-wave velocity                           [ m/s ]\n\
rho       density                                [ kg/m^3 ]\n\
\n\
--- end of file ---\n";

char *srcFormat =
"Format of source file\n\
\n\
--- beginning of file ---\n\
2              <- number of source points\n\
1              <- Number of pts in source template, must be 1\n\
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
Components can be Sx, Sz, Bulk and Sxz\n\
\n";


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



int main(int argc, char * const argv[])
{

    
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
    struct grid3d g;
    
    //
    // Source parameters
    //
    double f0=0.0;         // nominal frequency
    
	//
    // Motion variables
    //
	double *vx, *vy, *vz, *P;
    
	//
	// Material properties (input)
	//
	double *vp, *rho, *irhox, *irhoy, *irhoz;
	
    //
    // Memory variables
    //
    double *mem_dPdx, *mem_dPdy, *mem_dPdz;
    double *mem_dvdx, *mem_dvdy, *mem_dvdz;
    
    //
    //  CPML coefficients
    //
    double *ax,   *bx,   *kx,   *ay,   *by,   *ky,   *az,   *bz,   *kz;
    double *ax_h, *bx_h, *kx_h, *ay_h, *by_h, *ky_h, *az_h, *bz_h, *kz_h;
    
    
	params.ab = &g.ab;
	set_defaults_3d(&g, &params);
    process_args(argc, argv, &params);
    double dt = params.dt;
    strcpy(out.basename, params.basename);
    
	if ( verbose ) {
        fprintf(stdout, "\n *** a_iso_3d - Acoustic wave propagation in 3D isotropic media ***\n\n");
        if ( verbose > 1 ) {
            fprintf(stdout, "  Model file: \t%s\n", params.modelfile);
            fprintf(stdout, "  Source file:\t%s\n", params.sourcefile);
            fprintf(stdout, "  Output file:\t%s\n", params.outputfile);
            fprintf(stdout, "  Time window:\t%lg (s)\n", params.duration);
            fprintf(stdout, "  Time step:  \t%lg (ms)\n", params.dt * 1.e3);
        }
        fprintf(stdout, "Reading size of model ... ");
    }
    read_grid3d_params(params.modelfile, &g);
	
	if ( g.dx != g.dz ) {
		fprintf(stderr, "Error: grid cell size must be equal in X and Z, aborting.\n");
		abort();
	}
	
	const size_t np = g.ab.np;
	if ( verbose ) {
        fprintf(stdout, "done.\n");
        if ( verbose > 1 ) {
            fprintf(stdout, "  Model parameters: \t\tX\t\tY\t\tZ\n");
            fprintf(stdout, "        grid points \t\t%zd\t\t%zd\t\t%zd\n", g.nx, g.ny, g.nz);
            fprintf(stdout, "        grid spacing\t\t%lg\t\t%lg\t\t%lg\n", g.dx, g.dy, g.dz);
            fprintf(stdout, "        grid origin \t\t%lg\t\t%lg\t\t%lg\n", g.x0, g.y0, g.z0);
            fprintf(stdout, "        grid padding\t\t%zd\t\t%zd\t\t%zd\n", np, np, np);
        }
        fprintf(stdout, "Allocating memory ... ");
        fflush(stdout);
    }
	
    size_t nnn = g.nx2*g.ny2*g.nz2;
    
	vx =  (double *)malloc(nnn*sizeof(double)); // particle velocities
    vy =  (double *)malloc(nnn*sizeof(double));
    vz =  (double *)malloc(nnn*sizeof(double));
    P =   (double *)malloc(nnn*sizeof(double));
    vp =  (double *)malloc(nnn*sizeof(double));
    rho = (double *)malloc(nnn*sizeof(double));

    irhox = (double *)malloc(nnn*sizeof(double));
    irhoy = (double *)malloc(nnn*sizeof(double));
    irhoz = (double *)malloc(nnn*sizeof(double));
    
    mem_dPdx = (double *)malloc(nnn*sizeof(double));
    mem_dPdy = (double *)malloc(nnn*sizeof(double));
    mem_dPdz = (double *)malloc(nnn*sizeof(double));
    mem_dvdx = (double *)malloc(nnn*sizeof(double));
    mem_dvdy = (double *)malloc(nnn*sizeof(double));
    mem_dvdz = (double *)malloc(nnn*sizeof(double));
    
    ax = (double *)malloc((g.nx2)*sizeof(double));
    bx = (double *)malloc((g.nx2)*sizeof(double));
    kx = (double *)malloc((g.nx2)*sizeof(double));
    ay = (double *)malloc((g.ny2)*sizeof(double));
    by = (double *)malloc((g.ny2)*sizeof(double));
    ky = (double *)malloc((g.ny2)*sizeof(double));
    az = (double *)malloc((g.nz2)*sizeof(double));
    bz = (double *)malloc((g.nz2)*sizeof(double));
    kz = (double *)malloc((g.nz2)*sizeof(double));
    
    ax_h = (double *)malloc((g.nx2-1)*sizeof(double));
    bx_h = (double *)malloc((g.nx2-1)*sizeof(double));
    kx_h = (double *)malloc((g.nx2-1)*sizeof(double));
    ay_h = (double *)malloc((g.ny2-1)*sizeof(double));
    by_h = (double *)malloc((g.ny2-1)*sizeof(double));
    ky_h = (double *)malloc((g.ny2-1)*sizeof(double));
    az_h = (double *)malloc((g.nz2-1)*sizeof(double));
    bz_h = (double *)malloc((g.nz2-1)*sizeof(double));
    kz_h = (double *)malloc((g.nz2-1)*sizeof(double));
    
    
	//
    // Read in model
    //
    
    if ( verbose ) {
        fprintf(stdout, "done.\nReading properties of materials ... ");
        fflush(stdout);
    }
    read_model_a(params.modelfile, &g, vp, rho);
	
	double vpmax = 0.0;
    for (size_t i=0; i<g.nx; ++i) {
        for (size_t j=0; j<g.nz; ++j) {
			size_t ind = i*g.nz+j;
            vpmax = vp[ind]>vpmax ? vp[ind] : vpmax;
        }
    }
    
	//
	// averaging material properties at staggered nodes
	//
    for ( size_t i=0; i<g.nx2-1; ++i ) {
        for ( size_t j=0; j<g.ny2; ++j ) {
            for ( size_t k=0; k<g.nz2; ++k ) {
                size_t ind1 = (i*    g.ny2+j)*g.nz2+k;
                size_t ind2 = ((i+1)*g.ny2+j)*g.nz2+k;
                irhox[ind1] = dt * 2./(rho[ind1]+rho[ind2]);
            }
        }
    }
	for ( size_t i=0; i<g.nx2; ++i ) {
        for ( size_t j=0; j<g.ny2-1; ++j ) {
            for ( size_t k=0; k<g.nz2; ++k ) {
                size_t ind1 = (i*g.ny2+j  )*g.nz2+k;
                size_t ind2 = (i*g.ny2+j+1)*g.nz2+k;
                irhoy[ind1] = dt * 2./(rho[ind1]+rho[ind2]);
            }
        }
    }
	for ( size_t i=0; i<g.nx2; ++i ) {
        for ( size_t j=0; j<g.ny2; ++j ) {
            for ( size_t k=0; k<g.nz2-1; ++k ) {
                size_t ind1 = (i*g.ny2+j)*g.nz2+k;
                size_t ind2 = (i*g.ny2+j)*g.nz2+k+1;
                irhoz[ind1] = dt * 2./(rho[ind1]+rho[ind2]);
            }
        }
    }
    
    for ( size_t n=0; n<nnn; ++n ) {
        vp[n] = rho[n]*vp[n]*vp[n]*dt;
    }
    
	if ( dt*vpmax*sqrt(1./(g.dx*g.dx) + 1./(g.dy*g.dy) + 1./(g.dz*g.dz)) > 1. ) {
		
		fprintf(stderr,"\n *** Warning!  Time step too large\n");
		fprintf(stderr, "dt should be smaller that %lf\naborting ***\n",
				dt*vpmax*sqrt(1./(g.dx*g.dx) + 1./(g.dy*g.dy) + 1./(g.dz*g.dz)));
		abort();
	}
    //
    // Read in source params
    //
    if ( verbose ) {
        fprintf(stdout, "\nReading source parameters ... ");
        fflush(stdout);
    }
    read_source_3d(params.sourcefile, &g, &src);
	f0 = src.s[0].f;
	compute_src_fct(&src, params.dt, 1.0);
	
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
                fprintf(stdout, "    %zd - strength:   \t%lg MPa\n", ns+1, src.s[ns].A*fac);
                fprintf(stdout, "    %zd - coordinates:\t(%lg, %lg, %lg)\n", ns+1, src.s[ns].x, src.s[ns].y, src.s[ns].z);
            }
        }
        fprintf(stdout, "Reading output parameters ... ");
        fflush(stdout);
    }
    //
    // Read in output params
    //
    read_output_3d(params.outputfile, &g, &out);
    check_dt_trc(&out, dt);
    
    if ( params.segy == 1 ) {
        for ( size_t n=0; n<out.nrec; ++n ) {
            if ( out.r[n].type == TRACE ) {
                size_t nsamples = round( params.duration / out.r[n].dt );
                if ( NULL == ( out.r[n].data = (float *) malloc( nsamples * sizeof(float) ))) {
                    fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            }
        }
    }
    
    if ( verbose ) {
        fprintf(stdout, "done.\n");
        if ( verbose > 1 ) {
            fprintf(stdout, "There will be %zd record(s) on output\n", out.nrec);
            if ( params.segy == 1 ) {
                fprintf(stdout, "  Traces will be saved in SEG Y format in segy/%s.segy\n", params.basename);
            }
            fprintf(stdout, "  Snapshots will ");
            if ( !params.plotStrips ) fprintf(stdout, "not ");
            fprintf(stdout, "include absorbing strips.\n");
            for ( size_t nr=0; nr<out.nrec; ++nr ) {
                if ( out.r[nr].type == TRACE ) {
                    fprintf(stdout, "    %zd - Trace of %s at (%lg, %lg, %lg), ",
                            nr+1, component[out.r[nr].comp], out.r[nr].x, out.r[nr].y, out.r[nr].z);
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
    
	for (size_t n=0; n<nnn; ++n) {
        P[n] = vx[n] = vy[n] = vz[n] = 0.0;
        mem_dPdx[n] = mem_dPdy[n] = mem_dPdz[n] = 0.0;
        mem_dvdx[n] = mem_dvdy[n] = mem_dvdz[n] = 0.0;
    }
    
    for ( size_t n=0; n<g.nx2; ++n ) {
        ax[n] = 0.0;
        bx[n] = kx[n] = 1.0;
    }
    for ( size_t n=0; n<g.nx2-1; ++n ) {
        ax_h[n] = 0.0;
        bx_h[n] = kx_h[n] = 1.0;
    }
    for ( size_t n=0; n<g.ny2; ++n ) {
        ay[n] = 0.0;
        by[n] = ky[n] = 1.0;
    }
    for ( size_t n=0; n<g.ny2-1; ++n ) {
        ay_h[n] = 0.0;
        by_h[n] = ky_h[n] = 1.0;
    }
    for ( size_t n=0; n<g.nz2; ++n ) {
        az[n] = 0.0;
        bz[n] = kz[n] = 1.0;
    }
    for ( size_t n=0; n<g.nz2-1; ++n ) {
        az_h[n] = 0.0;
        bz_h[n] = kz_h[n] = 1.0;
    }
    
    const double pi = 4.0*atan(1.0);
    double alpha_max = pi*src.s[0].f;
    if ( g.ab.alpha_max == 0 ) alpha_max = 0.0;
    double da = alpha_max/np;
    double Lx = np*g.dx;
    double Ly = np*g.dy;
    double Lz = np*g.dz;
    double d0x = -(g.ab.pmlOrder + 1.) * g.ab.vmax * log(g.ab.Rc)/(2.*Lx);
    double d0y = -(g.ab.pmlOrder + 1.) * g.ab.vmax * log(g.ab.Rc)/(2.*Ly);
    double d0z = -(g.ab.pmlOrder + 1.) * g.ab.vmax * log(g.ab.Rc)/(2.*Lz);
    
    for ( size_t n=0; n<np; ++n ) {
        
        kx[n]   = kx[g.nx2-1-n]   = 1.+(g.ab.kappa_max-1.) * pow(((g.ab.np-n    )*g.dx/Lx), g.ab.pmlOrder);
        kx_h[n] = kx_h[g.nx2-2-n] = 1.+(g.ab.kappa_max-1.) * pow(((g.ab.np-n-0.5)*g.dx/Lx), g.ab.pmlOrder);
        ky[n]   = ky[g.ny2-1-n]   = 1.+(g.ab.kappa_max-1.) * pow(((g.ab.np-n    )*g.dy/Ly), g.ab.pmlOrder);
        ky_h[n] = ky_h[g.ny2-2-n] = 1.+(g.ab.kappa_max-1.) * pow(((g.ab.np-n-0.5)*g.dy/Ly), g.ab.pmlOrder);
        kz[n]   = kz[g.nz2-1-n]   = 1.+(g.ab.kappa_max-1.) * pow(((g.ab.np-n    )*g.dz/Lz), g.ab.pmlOrder);
        kz_h[n] = kz_h[g.nz2-2-n] = 1.+(g.ab.kappa_max-1.) * pow(((g.ab.np-n-0.5)*g.dz/Lz), g.ab.pmlOrder);
        
        double alpha = n*da;
        double alpha_h = (n+0.5)*da;
        
        double dx   = d0x * pow(((np-n    )*g.dx/Lx), g.ab.pmlOrder);
        double dx_h = d0x * pow(((np-n-0.5)*g.dx/Lx), g.ab.pmlOrder);

        bx[n]   = bx[g.nx2-1-n]   = exp( -(dx/kx[n]+alpha)*dt );
        bx_h[n] = bx_h[g.nx2-2-n] = exp( -(dx_h/kx_h[n]+alpha_h)*dt );
        
        ax[n]   = ax[g.nx2-1-n]   = dx  *(bx[n]  -1.)/(kx[n]  *(dx  +kx[n]  *alpha  ));
        ax_h[n] = ax_h[g.nx2-2-n] = dx_h*(bx_h[n]-1.)/(kx_h[n]*(dx_h+kx_h[n]*alpha_h));
        
        double dy   = d0y * pow(((np-n    )*g.dy/Ly), g.ab.pmlOrder);
        double dy_h = d0y * pow(((np-n-0.5)*g.dy/Ly), g.ab.pmlOrder);
        
        by[n]   = by[g.ny2-1-n]   = exp( -(dy/ky[n]+alpha)*dt );
        by_h[n] = by_h[g.ny2-2-n] = exp( -(dy_h/ky_h[n]+alpha_h)*dt );
        
        ay[n]   = ay[g.ny2-1-n]   = dy  *(by[n]  -1.)/(ky[n]  *(dy  +ky[n]  *alpha  ));
        ay_h[n] = ay_h[g.ny2-2-n] = dy_h*(by_h[n]-1.)/(ky_h[n]*(dy_h+ky_h[n]*alpha_h));
        
        double dz   = d0z * pow(((np-n    )*g.dz/Lz), g.ab.pmlOrder);
        double dz_h = d0z * pow(((np-n-0.5)*g.dz/Lz), g.ab.pmlOrder);
        
        bz[n]   = bz[g.nz2-1-n]   = exp( -(dz/kz[n]+alpha)*dt );
        bz_h[n] = bz_h[g.nz2-2-n] = exp( -(dz_h/kz_h[n]+alpha_h)*dt );
        
        az[n]   = az[g.nz2-1-n]   = dz  *(bz[n]  -1.)/(kz[n]  *(dz  +kz[n]  *alpha  ));
        az_h[n] = az_h[g.nz2-2-n] = dz_h*(bz_h[n]-1.)/(kz_h[n]*(dz_h+kz_h[n]*alpha_h));
       
//        printf("k = %lf, d = %lf, a = %lf, b = %lf\n", kx[n], dx, ax[n], bx[n]);
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
    double *data2 = NULL;
	for (size_t it=0; it<nsteps; ++it) {
        
		double t = dt*(it+1);
		if ( (it+1)%50 == 0 && verbose ) {
            fprintf(stdout, "  Iteration %zd, t = %g s\n", it+1, t); fflush(stdout);
        }
#pragma omp parallel default(shared)
        {
#pragma omp sections
            {

#pragma omp section
                for ( size_t i=0; i<g.nx2-1; ++i ) {
                    for ( size_t j=0; j<g.ny2; ++j ) {
                        for ( size_t k=0; k<g.nz2; ++k ) {
                            size_t ind1 = ((i  )*g.ny2+j)*g.nz2+k;
                            size_t ind2 = ((i+1)*g.ny2+j)*g.nz2+k;
                            
                            double dPdx = (P[ind2]-P[ind1])/g.dx;
                            mem_dPdx[ind1] = bx_h[i]*mem_dPdx[ind1] + ax_h[i]*dPdx;
                            dPdx = dPdx/kx_h[i] + mem_dPdx[ind1];
                            
                            vx[ind1] -= irhox[ind1]*dPdx;

                        }
                    }
                }
                
#pragma omp section
                for ( size_t i=0; i<g.nx2; ++i ) {
                    for ( size_t j=0; j<g.ny2-1; ++j ) {
                        for ( size_t k=0; k<g.nz2; ++k ) {
                            size_t ind1 = (i*g.ny2+j  )*g.nz2+k;
                            size_t ind2 = (i*g.ny2+j+1)*g.nz2+k;
                            
                            double dPdy = (P[ind2]-P[ind1])/g.dy;
                            mem_dPdy[ind1] = by_h[j]*mem_dPdy[ind1] + ay_h[j]*dPdy;
                            dPdy = dPdy/ky_h[j] + mem_dPdy[ind1];

                            vy[ind1] -= irhoy[ind1]*dPdy;
                        }
                    }
                }
                
#pragma omp section
                for ( size_t i=0; i<g.nx2; ++i ) {
                    for ( size_t j=0; j<g.ny2; ++j ) {
                        for ( size_t k=0; k<g.nz2-1; ++k ) {
                            size_t ind1 = (i*g.ny2+j)*g.nz2+k;
                            size_t ind2 = (i*g.ny2+j)*g.nz2+k+1;
                            
                            double dPdz = (P[ind2]-P[ind1])/g.dz;
                            mem_dPdz[ind1] = bz_h[k]*mem_dPdz[ind1] + az_h[k]*dPdz;
                            dPdz = dPdz/kz_h[k] + mem_dPdz[ind1];

                            vz[ind1] -= irhoz[ind1]*dPdz;
                        }
                    }
                }
            }  /*-- End of sections block --*/
			
			
#pragma omp sections
            {

#pragma omp section
                for ( size_t i=1; i<g.nx2; ++i ) {
                    for ( size_t j=1; j<g.ny2; ++j ) {
                        for ( size_t k=1; k<g.nz2; ++k ) {
                            size_t ind  = ((i  )*g.ny2+j  )*g.nz2+k;
                            size_t indx = ((i-1)*g.ny2+j  )*g.nz2+k;
                            size_t indy = ((i  )*g.ny2+j-1)*g.nz2+k;
                            size_t indz = ((i  )*g.ny2+j  )*g.nz2+k-1;
                            
                            double dvdx = (vx[ind]-vx[indx])/g.dx;
                            double dvdy = (vy[ind]-vy[indy])/g.dy;
                            double dvdz = (vz[ind]-vz[indz])/g.dz;
                            
                            mem_dvdx[ind] = bx[i]*mem_dvdx[ind] + ax[i]*dvdx;
                            dvdx = dvdx/kx[i] + mem_dvdx[ind];
                            
                            mem_dvdy[ind] = by[j]*mem_dvdy[ind] + ay[j]*dvdy;
                            dvdy = dvdy/ky[j] + mem_dvdy[ind];
                            
                            mem_dvdz[ind] = bz[k]*mem_dvdz[ind] + az[k]*dvdz;
                            dvdz = dvdz/kz[k] + mem_dvdz[ind];
                            
                            P[ind] -= vp[ind] * (dvdx + dvdy + dvdz);
                        }
                    }
                }
            }  /*-- End of sections block --*/
            
        } /*-- End of parallel region --*/
        
        // add source
        for (size_t ns=0; ns<src.nsrc; ++ns) {
            if ( it < src.s[ns].length ) {
                P[ src.s[ns].i ]    += dt * src.s[ns].fct[ it ] * 1.e6; // now in Pa
            }
        }
        
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
                case VY:
                    data = vy;
                    break;
                case VZ:
                    data = vz;
                    break;
                default:
                    break;
            }
			
            if ( ( fabs(dt2) < 0.99*dt ) ||
                ( ( t > out.r[nr].t0 && out.r[nr].dt > 0.0 && ( fabs(rem) < 1.e-8 ) ) ) ) {
				
                if ( out.r[nr].type == TRACE && params.segy==1 )
                    fill_data_segy(data, t, &out, nr);
                else if ( out.r[nr].type == TRACE && params.segy==0 )
                    write_trace(data, t, &out, nr);
                else if ( out.r[nr].type == SNAPSHOT ) {
                    
                    size_t iy = (src.s[0].y - g.y0)/g.dy;  // snapshot in the plane where src is
                    if ( data2 == NULL ) {
                        data2 = malloc(g.nx2*g.nz2*sizeof(double));
                    }
                    for ( size_t i=0, n=0; i<g.nx2; ++i ) {
                        for ( size_t k=0; k<g.nz2; ++k, ++n ) {
                            data2[n] = data[(i*g.ny2+iy)*g.nz2+k];
                        }
                    }
                    
					if ( verbose > 1 )
                        printf("Writing snapshot of %s, t = %g s\n", component[out.r[nr].comp], t);
                    write_snapshot3d_nc(data2, t, &g, &out, nr, params.plotStrips);
				}
            }
		}
	}
    
    if ( params.segy == 1 ) {
        save_segy(&params, &out, &src);
        for ( size_t n=0; n<out.nrec; ++n ) {
            if ( out.r[n].type == TRACE ) {
                free( out.r[n].data );
            }
        }
    }
    else
        for ( size_t n=0; n<out.nrec; ++n )
            if ( out.r[n].type == TRACE )
                fclose( out.r[n].fid );
    
	for (size_t n=0; n<src.nsrc; ++n) free( src.s[n].fct );
	free ( src.s );

	
	free( vx );
    free( vy );
    free( vz );
    free( P );
    
    free( vp );
    free( rho );
    free( irhox );
    free( irhoy );
    free( irhoz );
    
    free( mem_dPdx );
    free( mem_dPdy );
    free( mem_dPdz );
    free( mem_dvdx );
    free( mem_dvdy );
    free( mem_dvdz );
    
    free( ax );
    free( bx );
    free( kx );
    free( ay );
    free( by );
    free( ky );
    free( az );
    free( bz );
    free( kz );

    free( ax_h );
    free( bx_h );
    free( kx_h );
    free( ay_h );
    free( by_h );
    free( ky_h );
    free( az_h );
    free( bz_h );
    free( kz_h );
    
    free( data2 );

    if ( verbose ) fprintf(stdout, "\nEnd of computation\n");    
    return 0;

}

