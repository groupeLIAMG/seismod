//
//  e_cyl_iso.c
//
//  Elastic wave propagation in 2.5D isotropic media, on a staggered grid in
//     cylindrical coordinates with vertical symmetry axis
//
//  Created by Bernard Giroux on 15-10-14.
//
//

/*
 * Copyright (c) 2015, Bernard Giroux
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
 
 Reference papers:
 
 @Article{randall91,
 Title   = {Multipole borehole acoustic waveforms: Synthetic logs with beds and borehole washouts},
 Author  = {C. J. Randall and D. J. Scheibner and P. T. Wu},
 Journal = {Geophysics},
 Year    = {1991},
 Number  = {11},
 Pages   = {1757--1769},
 Volume  = {56},
 Doi     = {10.1190/1.1442988},
 Url     = {http://dx.doi.org/10.1190/1.1442988}
 }

 @Article{kurkjian86b,
 Title                    = {Acoustic multipole sources in fluid-filled boreholes},
 Author                   = {Andrew L. Kurkjian and Shu-Kong Chang},
 Journal                  = {Geophysics},
 Year                     = {1986},
 Number                   = {1},
 Pages                    = {148-163},
 Volume                   = {51},
 Doi                      = {10.1190/1.1442028},
 Url                      = {http://dx.doi.org/10.1190/1.1442028}
 }

 
 */

/*
 ------------------------------------------------------------------------------
 
 Usage
 
 e_cyl_iso [options] -p parameter_file.dat
 
 options are:
 -h Print this message
 -v verbose mode (type twice for increased verbosity)
 
 
 Discretization
 
 
 O---------#---------O
 |                   |
 |                   |
 |                   |
 |                   |
 *         @         *
 |                   |
 |                   |
 |                   |
 |                   |
 O---------#---------O

 
nodes         indices           variables
-----         -------           ---------
  O             i,j             tau_rz
  #           1+1/2,j           v_z, tau_thetaz
  *           i,j+1/2           v_r, tau_rtheta
  @         i+1/2,j+1/2         v_theta, tau_rr,tau_thetatheta, tau_zz
 
 At r=0, v_r and tau_rz are nil due to symmetry axis
 
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
1              # abs, apply CPML absorbing boundary (0 or 1, 0 by default)\n\
20             # nl, number of absorbing layers (20 by default)\n\
1              # plstr, plot absorbing strips in snapshots (0 or 1, 0 by default)\n\
0              # azimuthal_mode, azimuthal mode number (0 by default)\n\
--- end of file ---\n";

char *modelFormat =
"Format of model file\n\
\n\
--- beginning of file ---\n\
2 2                     <-  nr and nz, number of cells in r & z\n\
1.0 1.0                 <-  dr and dz, grid step size in r & z, in meters\n\
0.0 0.0                 <-  origin of the grid, in meters\n\
isotropic_elastic_cyl   <-  keyword for type of model\n\
V_p(1,1) V_s(1,1) rho(1,1)   <- properties for cell (1,1)\n\
V_p(2,1) V_s(2,1) rho(2,1)   <- properties for cell (2,1)\n\
...\n\
Vp(nr,1) Vs(nr,1) rho(nr,1)\n\
V_p(1,2) V_s(1,2) rho(1,2)\n\
V_p(2,2) V_s(2,2) rho(2,2)\n\
...\n\
V_p(nr,nz) V_s(nr,nz) rho(nr,nz)\n\
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
1D models can be input by giving values for just a horizontal profile; they\n\
will be duplicated automatically.  For example:\n\
\n\
--- beginning of file ---\n\
2 2            <-  nr and nz, number of cells in r & z\n\
1.0 1.0        <-  dr and dz, grid step size in r & z, in meters\n\
0.0 0.0        <-  origin of the grid, in meters\n\
isotropic_elastic      <-  keyword for type of model (isotropic, anisotropic_vti)\n\
V_p(1,1) V_s(1,1) rho(1,1)\n\
V_p(2,1) V_s(2,1) rho(2,1)\n\
...\n\
V_p(nr,1) V_s(nr,1) rho(nr,1)\n\
--- end of file ---\n";

char *srcFormat =
"Format of source file\n\
\n\
--- beginning of file ---\n\
1              <- number of source points\n\
\n\
Kurkjian       <- Source component of 1st source point\n\
1.0e-6         <- Peak volume change im m^3 at 1st source point\n\
50.0           <- Frequency in Hz at 1st source point\n\
50.0 45.0      <- r z coordinates at 1st source point in meters\n\
--- end of file ---\n\
\n\
Components can be Kurkjian, Fr, Fz, Ftheta\n\
\n\
}\n";


char *rcvFormat =
"Format of file defining the outputs\n\
\n\
--- beginning of file ---\n\
2              <- number of records\n\
\n\
Vz             <- particle velocity component of 1st record\n\
Trace          <- type of 1st record\n\
10.0 20.0      <- R Z coordinates at 1st record point in meters\n\
0.0            <- start time of 1st recording, in seconds\n\
2.0            <- time sampling of 1st record, in miliseconds\n\
\n\
Vr             <- particle velocity component of 2nd record\n\
Snapshot       <- type of 2nd record\n\
0.0            <- start time of snapshot(s) (2nd recording), in seconds\n\
2.0            <- time sampling of snapshot(s) 2nd record, in miliseconds.\n\
Give 0.0 or less for just one snapshot at time given above.\n\
--- end of file ---\n\
\n\
Admissible components are Vr, Vtheta, Vz, trr, trtheta, trz, tthetatheta, tthetaz, tzz, div, curl\n\
(div and curl are divergence and curl of particle velocity)\n\
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
Velocity components in traces are _not_ interpolated at exact R and Z\n\
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
    
    const char *typeSrc[] = { "Sx", "Sy", "Sz", "Sxy", "Sxz", "Syz", "Bulk", "Sf", "Bulk_s", "Fr", "Ftheta", "Fz", "Kurkjian" };
    const char *component[] = { "Vx", "Vy", "Vz", "qx", "qz", "tau_xx", "tau_zz", "tau_xz", "p", "tau_xy", "tau_yz", "divergence", "curl",
        "Vr", "Vtheta", "tau_rr", "tau_rtheta", "tau_rz", "tau_thetatheta", "tau_thetaz" };
    
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
    // Motion variables + material properties
    //
    struct variables_cyl v;
    
    //
    // Material properties (input)
    //
    double *vp, *vs, *rho, *mu_tmp;
    
    //
    // pml variables
    //
    
    params.ab = &g.ab;
    set_defaults(&g, &params);
    process_args(argc, argv, &params);
    double dt = params.dt;
    strcpy(out.basename, params.basename);
    
    if ( verbose ) {
        fprintf(stdout, "\n *** e_cyl_iso - Elastic wave propagation in 2.5D isotropic media in cylindrical coordinates ***\n\n");
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
            fprintf(stdout, "  Model parameters: \t\tR\t\tZ\n");
            fprintf(stdout, "        grid points \t\t%zd\t\t%zd\n", g.nx, g.nz);
            fprintf(stdout, "        grid spacing\t\t%lg\t\t%lg\n", g.dx, g.dz);
            fprintf(stdout, "        grid origin \t\t%lg\t\t%lg\n", g.x0, g.z0);
            fprintf(stdout, "        grid padding\t\t%zd\t\t%zd\n", g.ab.np, g.ab.np);
        }
        fprintf(stdout, "Allocating memory ... ");
        fflush(stdout);
    }
    
    size_t nnodes = g.nx2*g.nz2;
    if ( NULL == ( v.vr =  (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.vt =  (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.vz =  (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.trr = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.ttt = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.tzz = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.trz = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.trt = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.ttz = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.lij = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.l2mij = (double *)malloc(nnodes*sizeof(double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.mu =  (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.mui = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.muj = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.bi =  (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.bj =  (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( v.bij = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
    // temporary holders...
    mu_tmp = v.trr;
    vp = v.vr;
    vs = v.vz;
    rho = v.vt;
    
    
    //
    // Read in model
    //
    
    if ( verbose ) {
        fprintf(stdout, "done.\nReading properties of materials ... ");
        fflush(stdout);
    }
    read_model_e_cyl(params.modelfile, &g, vp, vs, rho);
    
    // properties are defined at i+1/2, j+1/2 (Randall et al, 1991)
    
    double vpmax = 0.0;
    for (size_t i=0; i<g.nx2; ++i) {
        for (size_t j=0; j<g.nz2; ++j) {
            size_t ind = i*g.nz2+j;
            vpmax = vp[ind]>vpmax ? vp[ind] : vpmax;
            v.lij[ind] = rho[ind]*(vp[ind]*vp[ind]-2.*vs[ind]*vs[ind]);
            mu_tmp[ind] = rho[ind]*vs[ind]*vs[ind];
            v.l2mij[ind] = v.lij[ind] + 2.0*mu_tmp[ind];
            v.bij[ind] = 1./rho[ind];
        }
    }
    
    //
    // averaging material properties
    //
    for (size_t i=0; i<g.nx2; ++i) {
        v.muj[i*g.nz2] = mu_tmp[i*g.nz2];
        v.bj[i*g.nz2] = v.bij[i*g.nz2];
        for (size_t j=1; j<g.nz2; ++j) {
            if (mu_tmp[i*g.nz2+j]==0.0 && mu_tmp[i*g.nz2+j-1]==0.0)
                v.muj[i*g.nz2+j] = 0.0;
            else
                v.muj[i*g.nz2+j] = 2.*mu_tmp[i*g.nz2+j]*mu_tmp[i*g.nz2+j-1]/(mu_tmp[i*g.nz2+j]+mu_tmp[i*g.nz2+j-1]);  // harmonic average
            double rho_tmp = (rho[i*g.nz2+j] + rho[i*g.nz2+j-1])/2.0; // arithmetic average
            v.bj[i*g.nz2+j] = 1./rho_tmp;
        }
    }
    for (size_t j=0; j<g.nz2; ++j) {
        v.mui[j] = mu_tmp[j];
        v.bi[j] = 1./rho[j];
        for (size_t i=1; i<g.nx2; ++i) {
            if (mu_tmp[i*g.nz2+j]==0.0 && mu_tmp[(i-1)*g.nz2+j]==0.0)
                v.mui[i*g.nz2+j] = 0.0;
            else
                v.mui[i*g.nz2+j] = 2.*mu_tmp[i*g.nz2+j]*mu_tmp[(i-1)*g.nz2+j]/(mu_tmp[i*g.nz2+j]+mu_tmp[(i-1)*g.nz2+j]);
            double rho_tmp = ( rho[i*g.nz2+j] + rho[(i-1)*g.nz2+j] )/2.0;
            v.bi[i*g.nz2+j] = 1./rho_tmp;
        }
    }
    v.mu[0] = mu_tmp[0];
    for (size_t j=1; j<g.nz2; ++j) {
        if (mu_tmp[j]==0.0 && mu_tmp[j-1]==0.0)
            v.mu[j] = 0.0;
        else
            v.mu[j] = 2.*mu_tmp[j]*mu_tmp[j-1]/(mu_tmp[j]+mu_tmp[j-1]);
    }
    for (size_t i=1; i<g.nx2; ++i) {
        
        if (mu_tmp[i*g.nz2]==0.0 && mu_tmp[(i-1)*g.nz2]==0.0)
            v.mu[i*g.nz2] = 0.0;
        else
            v.mu[i*g.nz2] = 2.*mu_tmp[i*g.nz2]*mu_tmp[(i-1)*g.nz2]/(mu_tmp[i*g.nz2]+mu_tmp[(i-1)*g.nz2]);
        
        for (size_t j=1; j<g.nz2; ++j) {
            
            if ( mu_tmp[i*g.nz2+j-1]*mu_tmp[(i-1)*g.nz2+j]*mu_tmp[(i-1)*g.nz2+j-1]==0.0 &&
                mu_tmp[i*g.nz2+j]*mu_tmp[(i-1)*g.nz2+j]*mu_tmp[(i-1)*g.nz2+j-1]==00 &&
                mu_tmp[i*g.nz2+j]*mu_tmp[i*g.nz2+j-1]*mu_tmp[(i-1)*g.nz2+j-1]==0.0 &&
                mu_tmp[i*g.nz2+j]*mu_tmp[i*g.nz2+j-1]*mu_tmp[(i-1)*g.nz2+j]==0.0 ) {
                v.mu[i*g.nz2+j] = 0.0;
            }
            else {
                v.mu[i*g.nz2+j] = 4.*mu_tmp[i*g.nz2+j]*mu_tmp[i*g.nz2+j-1]*mu_tmp[(i-1)*g.nz2+j]*mu_tmp[(i-1)*g.nz2+j-1] /
                (                  mu_tmp[i*g.nz2+j-1]*mu_tmp[(i-1)*g.nz2+j]*mu_tmp[(i-1)*g.nz2+j-1] +
                 mu_tmp[i*g.nz2+j]*                    mu_tmp[(i-1)*g.nz2+j]*mu_tmp[(i-1)*g.nz2+j-1] +
                 mu_tmp[i*g.nz2+j]*mu_tmp[i*g.nz2+j-1]*                      mu_tmp[(i-1)*g.nz2+j-1] +
                 mu_tmp[i*g.nz2+j]*mu_tmp[i*g.nz2+j-1]*mu_tmp[(i-1)*g.nz2+j]);
            }
        }
    }
    
    
//    write_field_nc(vp, "vp", "m/s", &g, &out, 1);
//    write_field_nc(vs, "vs", "m/s", &g, &out, 1);
//    write_field_nc(rho, "rho", "kg/m^3", &g, &out, 1);
//    
//    write_field_nc(v.bi,  "bi", "m^3/kg", &g, &out, 1);
//    write_field_nc(v.bj,  "bj", "m^3/kg", &g, &out, 1);
//    write_field_nc(v.bij, "bij", "m^3/kg", &g, &out, 1);
//    write_field_nc(v.mu,  "mu", "Pa", &g, &out, 1);
//    write_field_nc(v.mui, "mui", "Pa", &g, &out, 1);
//    write_field_nc(v.muj, "muj", "Pa", &g, &out, 1);
//    write_field_nc(v.lij, "lij", "Pa", &g, &out, 1);
//    write_field_nc(v.l2mij, "l2mij", "Pa", &g, &out, 1);

    
    
    for ( size_t n=0; n<nnodes; ++n ) {
        
        if ( isfinite(v.mu[n])==0 ) {
            fprintf(stderr, "Error, mu not a finite value at node %zd.  Aborting", n);
            abort();
        }
        if ( isfinite(v.mui[n])==0 ) {
            fprintf(stderr, "Error, mu not a finite value at node %zd (staggered i+1/2).  Aborting", n);
            abort();
        }
        if ( isfinite(v.muj[n])==0 ) {
            fprintf(stderr, "Error, mu not a finite value at node %zd (staggered j+1/2).  Aborting", n);
            abort();
        }
        if ( isfinite(v.bi[n])==0 ) {
            fprintf(stderr, "Error, 1/rho not a finite value at node %zd (staggered i+1/2).  Aborting", n);
            abort();
        }
        if ( isfinite(v.bj[n])==0 ) {
            fprintf(stderr, "Error, 1/rho not a finite value at node %zd (staggered j+1/2).  Aborting", n);
            abort();
        }
        if ( isfinite(v.bij[n])==0 ) {
            fprintf(stderr, "Error, 1/rho not a finite value at node %zd (staggered i+1/2,j+1/2).  Aborting", n);
            abort();
        }
        if ( isfinite(v.lij[n])==0 ) {
            fprintf(stderr, "Error, Lamé lambda not a finite value at node %zd (staggered i+1/2,j+1/2).  Aborting", n);
            abort();
        }
        if ( isfinite(v.l2mij[n])==0 ) {
            fprintf(stderr, "Error, lambda+2mu not a finite value at node %zd (staggered i+1/2,j+1/2).  Aborting", n);
            abort();
        }
    }
    
    double Delta = 0.5*(g.dx+g.dz);
    if ( dt > Delta / ( sqrt(2.0) * vpmax ) ) {
        fprintf(stderr, "\n *** Warning!  Time step too large, aborting ***\n");
        abort();
    }

    //
    // Read in source params
    //
    if ( verbose ) {
        fprintf(stdout, "done.\nReading source parameters ... ");
        fflush(stdout);
    }
    read_source(params.sourcefile, &g, &src);
    f0 = src.s[0].f;
    if ( src.s[0].type == KURKJIAN ) {
        compute_kurkjian(&src, &v, &g, params.dt, params.n);
    } else {
        compute_src_fct(&src, params.dt, 1.0);
    }
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
                if ( src.s[ns].type == KURKJIAN )
                    fprintf(stdout, "    %zd - Peak vol change:\t%lg m^3\n", ns+1, src.s[ns].A);
                else
                    fprintf(stdout, "    %zd - strength:   \t%lg MPa\n", ns+1, src.s[ns].A*1.e3*fac);
                fprintf(stdout, "    %zd - coordinates:\t(%lg, %lg)\n", ns+1, src.s[ns].x, src.s[ns].z);
            }
        }
    }
    
    //
    // CPML stuff
    //
    if ( verbose ) fprintf(stdout, "done.\nApplying PML absorbing boundaries (%zd layers) ", g.ab.np);
    struct mem_cpml_cyl mem;
    struct fac_cpml_cyl fp;
    alloc_cpml_cyl(&fp, &mem, &g, params.n);
    const double pi = 4.0*atan(1.0);
    double alpha_max = pi*src.s[0].f;
    if ( g.ab.alpha_max == 0 ) alpha_max = 0.0;
    compute_cpml_cyl(&fp, &g, params.dt, alpha_max, params.iwipe);
    
    //
    // Read in output params
    //
    if ( verbose ) {
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
    double *div=NULL;
    double *curl=NULL;
    for ( size_t nr=0; nr<out.nrec; ++nr ) {
        if (out.r[nr].comp == DIV) {
            if ( NULL == ( div = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            for (size_t n=0; n<nnodes; ++n) { div[n] = 0.0; }
        }
        else if (out.r[nr].comp == CURL) {
            if ( NULL == ( curl = (double *)malloc(nnodes*sizeof(double))))   { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            for (size_t n=0; n<nnodes; ++n) { curl[n] = 0.0; }
        }
    }
    
    for (size_t n=0; n<nnodes; ++n) {
        v.vr[n] = v.vt[n] = v.vz[n] = v.trr[n] = v.ttt[n] = v.tzz[n] = v.trz[n] = v.trt[n] = v.ttz[n] = 0.0;
    }
    
    
    nsteps = params.duration / dt;
    if ( nsteps*dt < params.duration ) nsteps++;
    
#ifdef _OPENMP
    if ( verbose ) {
        fprintf(stdout, "Will be running threaded (OpenMP) version with 4 threads\n");
    }
    omp_set_num_threads( 4 );
#endif
    
    if ( verbose ) {
        fprintf(stdout, "Starting main loop (%zd iterations)\n", nsteps);
        fflush(stdout);
    }
    if ( params.n > 0 ) {
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
                    // Vr
                    // ------------------------------------------------------------
#pragma omp section
                    update_vr_cyl_n(&v, &g, &mem, &fp, dt, params.n);
                    
                    // ------------------------------------------------------------
                    // Vtheta
                    // ------------------------------------------------------------
#pragma omp section
                    update_vt_cyl_n(&v, &g, &mem, &fp, dt, params.n);
                    
                    // ------------------------------------------------------------
                    // Vz
                    // ------------------------------------------------------------
#pragma omp section
                    update_vz_cyl_n(&v, &g, &mem, &fp, dt, params.n);
                    
                } /*-- End of sections block --*/
                
                add_force_src_cyl(&src, &v, dt, it);
                
#pragma omp sections
                {
                    // ------------------------------------------------------------
                    // tau_rr, tau_\alpha\alpha and tau_zz
                    // ------------------------------------------------------------
#pragma omp section
                    update_taa_cyl_n(&v, &g, &mem, &fp, dt, params.n);
                    
                    // ------------------------------------------------------------
                    // tau_rz
                    // ------------------------------------------------------------
#pragma omp section
                    update_trz_cyl(&v, &g, &mem, &fp, dt);
                    
                    // ------------------------------------------------------------
                    // tau_r\theta
                    // ------------------------------------------------------------
#pragma omp section
                    update_trt_cyl_n(&v, &g, &mem, &fp, dt, params.n);
                    
                    // ------------------------------------------------------------
                    // tau_theta z
                    // ------------------------------------------------------------
#pragma omp section
                    update_ttz_cyl_n(&v, &g, &mem, &fp, dt, params.n);
                    
                } /*-- End of sections block --*/
            } /*-- End of parallel region --*/
            
            if ( src.s[0].type == KURKJIAN ) {
                add_kurkjian_src(&src, &v, dt, it);
            }
            
            // ------------------------------------------------------------
            // write records
            // ------------------------------------------------------------
            for ( size_t nr=0; nr<out.nrec; ++nr ) {
                
                double dt2 = t-out.r[nr].t0+1.e-15;
                double rem = fmod(1.e6*dt2, 1.e6*out.r[nr].dt);
                double *data=NULL;
                switch (out.r[nr].comp) {
                    case VR:
                        data = v.vr;
                        break;
                    case VT:
                        data = v.vt;
                        break;
                    case VZ:
                        data = v.vz;
                        break;
                    case TRR:
                        data = v.trr;
                        break;
                    case TRT:
                        data = v.trt;
                        break;
                    case TRZ:
                        data = v.trz;
                        break;
                    case TTT:
                        data = v.ttt;
                        break;
                    case TTZ:
                        data = v.ttz;
                        break;
                    case TZZ:
                        data = v.tzz;
                        break;
                    case DIV:
                        compute_div_cyl(div, &g, &v);
                        data = div;
                        break;
                    case CURL:
                        compute_curl_cyl(curl, &g, &v, params.n);
                        data = curl;
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
    } else {
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
                    // Vr
                    // ------------------------------------------------------------
#pragma omp section
                    update_vr_cyl(&v, &g, &mem, &fp, dt);
                    
                    // ------------------------------------------------------------
                    // Vtheta
                    // ------------------------------------------------------------
#pragma omp section
                    update_vt_cyl(&v, &g, &mem, &fp, dt);
                    
                    // ------------------------------------------------------------
                    // Vz
                    // ------------------------------------------------------------
#pragma omp section
                    update_vz_cyl(&v, &g, &mem, &fp, dt);
                    
                } /*-- End of sections block --*/
                
                add_force_src_cyl(&src, &v, dt, it);
                
#pragma omp sections
                {
                    // ------------------------------------------------------------
                    // tau_rr, tau_\alpha\alpha and tau_zz
                    // ------------------------------------------------------------
#pragma omp section
                    update_taa_cyl(&v, &g, &mem, &fp, dt);
                    
                    // ------------------------------------------------------------
                    // tau_rz
                    // ------------------------------------------------------------
#pragma omp section
                    update_trz_cyl(&v, &g, &mem, &fp, dt);
                    
                    // ------------------------------------------------------------
                    // tau_r\theta
                    // ------------------------------------------------------------
#pragma omp section
                    update_trt_cyl(&v, &g, &mem, &fp, dt);
                    
                    // ------------------------------------------------------------
                    // tau_theta z
                    // ------------------------------------------------------------
#pragma omp section
                    update_ttz_cyl(&v, &g, &mem, &fp, dt);
                    
                } /*-- End of sections block --*/
            } /*-- End of parallel region --*/
            
            if ( src.s[0].type == KURKJIAN ) {
                add_kurkjian_src(&src, &v, dt, it);
            }

            // ------------------------------------------------------------
            // write records
            // ------------------------------------------------------------
            for ( size_t nr=0; nr<out.nrec; ++nr ) {
                
                double dt2 = t-out.r[nr].t0+1.e-15;
                double rem = fmod(1.e6*dt2, 1.e6*out.r[nr].dt);
                double *data=NULL;
                switch (out.r[nr].comp) {
                    case VR:
                        data = v.vr;
                        break;
                    case VT:
                        data = v.vt;
                        break;
                    case VZ:
                        data = v.vz;
                        break;
                    case TRR:
                        data = v.trr;
                        break;
                    case TRT:
                        data = v.trt;
                        break;
                    case TRZ:
                        data = v.trz;
                        break;
                    case TTT:
                        data = v.ttt;
                        break;
                    case TTZ:
                        data = v.ttz;
                        break;
                    case TZZ:
                        data = v.tzz;
                        break;
                    case DIV:
                        compute_div_cyl(div, &g, &v);
                        data = div;
                        break;
                    case CURL:
                        compute_curl_cyl(curl, &g, &v, params.n);
                        data = curl;
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
    }
    
    free( v.vr );
    free( v.vt );
    free( v.vz );
    free( v.trr );
    free( v.ttt );
    free( v.tzz );
    free( v.trz );
    free( v.trt );
    free( v.ttz );
    free( v.lij );
    free( v.l2mij );
    free( v.mu );
    free( v.mui );
    free( v.muj );
    free( v.bi );
    free( v.bj );
    free( v.bij );
    
    free_cpml_cyl(&fp, &mem);
    
    if ( verbose ) fprintf(stdout, "\nEnd of computation\n");    
    return 0;
}
