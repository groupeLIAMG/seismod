//
//  ve_vti_pml.c
//
//  Created by Bernard Giroux on 2012-12-03.
//
//

//  ViscoElastic wave propagation in 2D VTI anisotropic media, for L
//  attenuation quasi-dilatational mechanism and 1 shear relaxation mechanism,
//  on a staggered grid
//
//  Code specs:  - language: ANSI C99
//               - compiled with intel compiler 13.0 on a mac running OSX 10.8
//               - external libraries: fftw ( http://www.fftw.org ),
//                       netCDF ( http://www.unidata.ucar.edu/software/netcdf/ )
//

/*
 *
 *   Reference papers
 *
 
 @ARTICLE{carcione95c,
 author = {Jose M. Carcione},
 title = {Constitutive model and wave equations for linear, viscoelastic, anisotropic
 media},
 journal = {Geophysics},
 year = {1995},
 volume = {60},
 pages = {537-548},
 number = {2},
 doi = {10.1190/1.1443791}
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
 doi = {10.1046/j.1365-2478.1998.00087.x}
 }

 @ARTICLE{carcione99c,
 author = {Jos\'{e} M. Carcione},
 title = {Staggered mesh for the anisotropic and viscoelastic wave equation},
 journal = {Geophysics},
 year = {1999},
 volume = {64},
 pages = {1863-1866},
 number = {6},
 doi = {10.1190/1.1444692}
 }

 ------------------------------------------------------------------------------
 
 Format of model file
 
 --- beginning of file ---
 2 2                  <-  nx and nz, number of nodes in x & z
 1.0 1.0              <-  dx and dz, grid step size in x & z, in meters
 0.0 0.0              <-  origin of the grid, in meters
 vti_viscoelastic     <-  keyword for type of model
 2                    <-  Number of relaxation mechanisms, quasi-dilatational mode
 K(1,1) mu(1,1) rho(1,1) epsilon(1,1) delta(1,1) Q1(1,1,1) f1(1,1,1) Q1(2,1,1) f1(2,1,1) Q2(1,1) f2(1,1)
 K(1,2) mu(1,2) rho(1,2) epsilon(1,2) delta(1,2) Q1(1,1,2) f1(1,1,2) Q1(2,1,2) f1(2,1,2) Q2(1,2) f2(1,2)
 K(2,1) mu(2,1) rho(2,1) epsilon(2,1) delta(2,1) Q1(1,2,1) f1(1,2,1) Q1(2,2,1) f1(2,2,1) Q2(2,1) f2(2,1)
 K(2,2) mu(2,2) rho(2,2) epsilon(2,2) delta(2,2) Q1(1,2,2) f1(1,2,2) Q1(2,2,2) f1(2,2,2) Q2(2,2) f2(2,2)
 --- end of file ---
 
 Variable  Description                                                     Units
 
 K         bulk modulus                                                  [ GPa ]
 mu        shear modulus of the matrix                                   [ GPa ]
 rho       density                                                    [ kg/m^3 ]
 epsilon   Thomsen anisotropy parameter                                    [ - ]
 delta     Thomsen anisotropy parameter                                    [ - ]
 Q1        Q for quasi-dilatational mode                                   [ - ]
 f1        Relaxation frequencies - Q1                                    [ Hz ]
 Q2        Q for shear mode                                                [ - ]
 f2        Relaxation frequency - Q2                                      [ Hz ]

 
 ------------------------------------------------------------------------------
 
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "io_utils.h"
#include "pml.h"
#include "propagate.h"
#include "structs.h"
#include "src.h"

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
0              # simulateCMP, for 1D model save SEGY as a CMP\n\
0              # check_model, compute phase vel and attenuation (0, 1 (display) or 2 (store), 0 by default)\n\
1              # abs, apply absorbing boundary (0 or 1, 0 by default)\n\
20             # nl, number of absorbing layers (20 by default)\n\
1              # plstr, plot absorbing strips in snapshots (0 or 1, 0 by default)\n\
7              # kappa_max, PML parameter (1 by default)\n\
1              # alpha_max, PML parameter.  If 0, alpha = 1 everywhere in the PML; if 1 alpha_max = pi*f0\n\
0.001          # Rc, theoretical coeff of reflection of PML layer (1e-5 by default)\n\
1              # saveEnergy, save E as fct of time\n\
0/500/0/500    # EnergyROI, region of interest for saving E\n\
1000           # chkpt_inc, checkpoint saving increment (0 by default, no checkpoint saving)\n\
--- end of file ---\n";

char *modelFormat =
"Format of model file\n\
\n\
--- beginning of file ---\n\
2 2                  <-  nx and nz, number of nodes in x & z\n\
1.0 1.0              <-  dx and dz, grid step size in x & z, in meters\n\
0.0 0.0              <-  origin of the grid, in meters\n\
vti_viscoelastic     <-  keyword for type of model\n\
2                    <-  Number of relaxation mechanisms, quasi-dilatational mode\n\
Vp0(1,1) Vs0(1,1) rho(1,1) epsilon(1,1) delta(1,1) Q1(1,1,1) f1(1,1,1) Q1(2,1,1) f1(2,1,1) Q2(1,1) f2(1,1)\n\
Vp0(1,2) Vs0(1,2) rho(1,2) epsilon(1,2) delta(1,2) Q1(1,1,2) f1(1,1,2) Q1(2,1,2) f1(2,1,2) Q2(1,2) f2(1,2)\n\
Vp0(2,1) Vs0(2,1) rho(2,1) epsilon(2,1) delta(2,1) Q1(1,2,1) f1(1,2,1) Q1(2,2,1) f1(2,2,1) Q2(2,1) f2(2,1)\n\
Vp0(2,2) Vs0(2,2) rho(2,2) epsilon(2,2) delta(2,2) Q1(1,2,2) f1(1,2,2) Q1(2,2,2) f1(2,2,2) Q2(2,2) f2(2,2)\n\
...\n\
K(nx,nz) mu(nx,nz) rho(nx,nz) epsilon(nx,nz) delta(nx,nz) Q1(1,nx,nz) f1(1,nx,nz) Q1(2,nx,nz) f1(2,nx,nz) Q2(nx,nz) f2(nx,nz)\n\
--- end of file ---\n\
\n\
Variable  Description                                                     Units\n\
\n\
Vp0       vertical P-wave velocity                                      [ m/s ]\n\
Vs0       vertical S-wave velocity                                      [ m/s ]\n\
rho       density                                                    [ kg/m^3 ]\n\
epsilon   Thomsen anisotropy parameter                                    [ - ]\n\
delta     Thomsen anisotropy parameter                                    [ - ]\n\
Q1        Q for quasi-dilatational mode                                   [ - ]\n\
f1        Relaxation frequencies - Q1                                    [ Hz ]\n\
Q2        Q for shear mode                                                [ - ]\n\
f2        Relaxation frequency - Q2                                      [ Hz ]\n\
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
vti_viscoelastic  <-  keyword for type of model\n\
2                 <-  Number of relaxation mechanisms, quasi-dilatational mode\n\
Vp0(1,1) Vs0(1,1) rho(1,1) epsilon(1,1) delta(1,1) Q1(1,1,1) f1(1,1,1) Q1(2,1,1) f1(2,1,1) Q2(1,1) f2(1,1)\n\
Vp0(1,2) Vs0(1,2) rho(1,2) epsilon(1,2) delta(1,2) Q1(1,1,2) f1(1,1,2) Q1(2,1,2) f1(2,1,2) Q2(1,2) f2(1,2)\n\
...\n\
Vp0(1,nz) Vs0(1,nz) rho(1,nz) epsilon(1,nz) delta(1,nz) Q1(1,1,nz) f1(1,1,nz) Q1(2,1,nz) f1(2,1,nz) Q2(1,nz) f2(1,nz)\n\
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
Components can be Sx, Sz, Bulk and Sxz\n\
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



int main(int argc, char * argv[])
{

    struct inputParams params;
    struct computationVariablesVE_VTI c;
    struct fac_pml fp1, fp23, fp4;
	struct mem_pml mem[3];
    struct materialPropertiesVE_VTI mp;
    struct fftw_data fd;
	struct sourceParams src;
    struct outputParams out;
//	struct saveEnergyVTI se;

	// -------------------------------------------------------------------------
    //
    // variables
    //
    // -------------------------------------------------------------------------
	
    const double pi = 4.0*atan(1.0);
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
    double *v_x;       // particle velocity along x
    double *v_z;       // particle velocity along z

    double *W;
    double *Ws;
    double *Wtmp;
    double *delta;
	
	size_t itstart = 0;
    
    params.ab = &g.ab;
	set_defaults(&g, &params);
    process_args(argc, argv, &params);

    strcpy(out.basename, params.basename);
    
    if ( verbose ) {
        fprintf(stdout, "\n *** ve_vti - ViscoElastic wave propagation in 2D VTI anisotropic media ***\n\n");
        if ( verbose > 1 ) {
            fprintf(stdout, "  Model file: \t%s\n", params.modelfile);
            fprintf(stdout, "  Source file:\t%s\n", params.sourcefile);
            fprintf(stdout, "  Output file:\t%s\n", params.outputfile);
            fprintf(stdout, "  Time window:\t%lg (s)\n", params.duration);
            fprintf(stdout, "  Time step:  \t%lg (ms)\n", params.dt * 1.e3);
//            if ( params.chkpt_inc > 0 )
//                fprintf(stdout, "  Checkpointing every %d iterations applied\n", params.chkpt_inc );
        }
        fprintf(stdout, "Reading size of model ... ");
    }
    read_grid_params_ve(params.modelfile, &g, &mp);
    
    if ( verbose ) {
        fprintf(stdout, "done.\n");
        if ( verbose > 1 ) {
            fprintf(stdout, "  Model parameters: \t\tX\t\tZ\n");
            fprintf(stdout, "        grid points \t\t%zd\t\t%zd\n", g.nx, g.nz);
            fprintf(stdout, "        grid spacing\t\t%lg\t\t%lg\n", g.dx, g.dz);
            fprintf(stdout, "        grid origin \t\t%lg\t\t%lg\n", g.x0, g.z0);
            fprintf(stdout, "        Nb of quasi-dilatational relaxation mechanisms\t%d\n", mp.L);
            fprintf(stdout, "        pml padding\t\t%zd\t\t%zd\n", g.ab.np, g.ab.np);
            fprintf(stdout, "        pml - Rc       \t\t%lg\n", g.ab.Rc);
            fprintf(stdout, "        pml - kappa_max\t\t%lg\n", g.ab.kappa_max);
        }
        fprintf(stdout, "Allocating memory ... ");
        fflush(stdout);
    }
	c.dt = params.dt;
	
	//
    //  Memory allocation
    //
    
    size_t nnodes = g.nx2 * g.nz2;
    
    if ( NULL == ( mp.K       = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.mu      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.rho     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.epsilon = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.delta   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.Q2      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.f2      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.Q1      = (double **) malloc( mp.L * sizeof(double*) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( mp.f1      = (double **) malloc( mp.L * sizeof(double*) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    for ( size_t n=0; n<mp.L; ++n ) {
        if ( NULL == ( mp.Q1[n] = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( mp.f1[n] = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    }

    if ( NULL == ( c.c11        = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.c13        = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.c33        = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.c55        = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.K0         = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.c55_0      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.c55_0_ij   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.epsilon2   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.epsilon3   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.tau_s_2    = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.tau_s_2_ij = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.eta2       = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.eta2_ij    = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.rho_i      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.rho_j      = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.epsilon1   = (double **) malloc( mp.L * sizeof(double*) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.tau_s_1    = (double **) malloc( mp.L * sizeof(double*) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( c.eta1       = (double **) malloc( mp.L * sizeof(double*) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    for ( size_t n=0; n<mp.L; ++n ) {
        if ( NULL == ( c.epsilon1[n] = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( c.tau_s_1[n]  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( c.eta1[n]     = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    }
    
    if ( NULL == ( W      = (double *) malloc( (7+mp.L)*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( Ws     = (double *) malloc( (7+mp.L)*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( Wtmp   = (double *) malloc( (7+mp.L)*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( delta  = (double *) malloc( (7+mp.L)*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	

    alloc_cpml_ve(&fp1, &fp23, &fp4, mem, &g);


    //
    // Read in model
    //
    
    if ( verbose ) {
        fprintf(stdout, "done.\nReading properties of materials ... ");
        fflush(stdout);
    }
    read_model_ve(params.modelfile, &g, &mp);
    
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
            fprintf(stdout, "  Source has %zd component(s)\n", src.nsrc);
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
    check_dt_trc(&out, params.dt);
    
    if ( params.segy == 1 ) {
        for ( size_t n=0; n<out.nrec; ++n ) {
            if ( out.r[n].type == TRACE ) {
                size_t nsamples = round( params.duration / out.r[n].dt );
                if ( NULL == ( out.r[n].data = (float *) malloc( nsamples * sizeof(float) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
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
    }
//    if ( params.saveEnergy == 1 ) {
//        struct stat sb;
//        if ( stat("E", &sb) < 0 ) {
//            if ( errno == ENOENT ) {
//                mkdir("E", S_IRWXU);
//            }
//        }
//        char fname[80];
//        sprintf(fname, "E/%s_E.dat", params.basename );
//        se.fid = fopen( fname, "w+" );
//        if ( verbose ) {
//            fprintf(stdout, "Saving kinetic energy as fct of time in %s\n", fname);
//            if (params.roiExmin<g.x0 || params.roiExmax>(g.x0+(g.nx-1)*g.dx) ||
//                params.roiEzmin<g.z0 || params.roiEzmax>(g.z0+(g.nz-1)*g.dz)) {
//                fprintf(stderr, "Region of interest for saving energy not fitting in the modeling grid\n");
//                abort();
//            }
//            
//            fprintf(stdout, "  Region of interest is %lg/%lg/%lg/%lg\n",
//                    params.roiExmin, params.roiExmax, params.roiEzmin, params.roiEzmax);
//        }
//        se.i1E = lround( (params.roiExmin-g.x0)/g.dx );
//        se.i2E = lround( (params.roiExmax-g.x0)/g.dx );
//        se.j1E = lround( (params.roiEzmin-g.z0)/g.dz );
//        se.j2E = lround( (params.roiEzmax-g.z0)/g.dz );
//    }    

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
        mp.K[n]   *= 1000.;          // now in MPa
        mp.mu[n]  *= 1000.;          // now in MPa
        mp.rho[n] *= 1.e-6;          // density are now in Gg/m^3 (Mega. kg/m^3)
    }
    
    //
    // Compute elastic tensor components & other physical params
    //
    for ( size_t n=0; n<nnodes; ++n ) {
        c.c33[n] = mp.K[n] + 4./3.*mp.mu[n];
        double c55 = mp.mu[n];
        c.c11[n] = 2. * mp.epsilon[n] * c.c33[n] + c.c33[n];
        c.c13[n] = sqrt( 2. * mp.delta[n] * c.c33[n] * (c.c33[n] - c55) +
                        (c.c33[n] - c55) * (c.c33[n] - c55) ) - c55;
        for ( size_t m=0; m<mp.L; ++m ) {
            double tau_e = (sqrt(mp.Q1[m][n]*mp.Q1[m][n]+1)+1.) / (2.*pi*mp.f1[m][n]*mp.Q1[m][n]);
            c.tau_s_1[m][n] = (sqrt(mp.Q1[m][n]*mp.Q1[m][n]+1)-1.) / (2.*pi*mp.f1[m][n]*mp.Q1[m][n]);
            c.eta1[m][n] = c.tau_s_1[m][n] / tau_e;
        }
        double tau_e = (sqrt(mp.Q2[n]*mp.Q2[n]+1)+1.) / (2.*pi*mp.f2[n]*mp.Q2[n]);
        c.tau_s_2[n] = (sqrt(mp.Q2[n]*mp.Q2[n]+1)-1.) / (2.*pi*mp.f2[n]*mp.Q2[n]);
        c.eta2[n] = c.tau_s_2[n] / tau_e;
    
        c.c55_0[n] = c55*c.eta2[n];
        double D = 0.5*(c.c11[n]+c.c33[n]);
        double K = D - c55;
        double eta1 = 0.0;
        for ( size_t m=0; m<mp.L; ++m ) {
            eta1 += c.eta1[m][n];
        }
        
        double c11_0 = c.c11[n] - D + K*eta1 + c55*c.eta2[n];
        double c33_0 = c.c33[n] - D + K*eta1 + c55*c.eta2[n];
//        double c13_0 = c.c13[n] + c55*(2.-c.eta2[n]);
        
        c.K0[n] = 0.5*(c11_0+c33_0) - c.c55_0[n];
    }
    
    //
    //  Averaging of parameters at staggered nodes, at (i+1/2,j)
    //
    for ( size_t i=0; i<g.nx2-1; ++i ) {
        for ( size_t j=0; j<g.nz2; ++j ) {
            
            size_t ij = i*g.nz2+j;
            size_t iij = (i+1)*g.nz2+j;
            
            c.rho_i[ij]   = 0.5 * ( mp.rho[ij] + mp.rho[iij] );
        }
    }

    //
    //  Averaging of parameters at staggered nodes, at (i,j+1/2)
    //
    for ( size_t i=0; i<g.nx2; ++i ) {
        for ( size_t j=0; j<g.nz2-1; ++j ) {
            
            size_t ij = i*g.nz2+j;
            size_t ijj = ij+1;
            
            c.rho_j[ij]   = 0.5 * ( mp.rho[ij] + mp.rho[ijj] );
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
            
            c.tau_s_2_ij[ij] = 0.25*(c.tau_s_2[ij] +
                                     c.tau_s_2[(i+1)*g.nz2+j] +
                                     c.tau_s_2[i*g.nz2+j+1] +
                                     c.tau_s_2[(i+1)*g.nz2+j+1]);
            
            c.eta2_ij[ij] = 0.25*(c.eta2[ij] +
                                  c.eta2[(i+1)*g.nz2+j] +
                                  c.eta2[i*g.nz2+j+1] +
                                  c.eta2[(i+1)*g.nz2+j+1]);
        }
    }

    //  fill end of grid
    for ( size_t i=g.nx2-1, j=0; j<g.nz2; ++j ) {
        
        size_t ij  = (i-1)*g.nz2+j;
        size_t iij = i*g.nz2+j;
        
        c.rho_i[iij]      = c.rho_i[ij];
        c.c55[iij]        = c.c55[ij];
        c.tau_s_2_ij[iij] = c.tau_s_2_ij[ij];
        c.eta2_ij[iij]    = c.eta2_ij[ij];
    }
    for ( size_t i=0, j=g.nz2-1; i<g.nx2; ++i ) {
        
        size_t ij  = i*g.nz2+j-1;
        size_t ijj = i*g.nz2+j;
        
        c.rho_j[ijj]      = c.rho_j[ij];
        c.c55[ijj]        = c.c55[ij];
        c.tau_s_2_ij[ijj] = c.tau_s_2_ij[ij];
        c.eta2_ij[ijj]    = c.eta2_ij[ij];        
    }
    
//    if ( params.saveEnergy == 1 ) {
//        size_t npts = (se.i2E-se.i1E+1)*(se.j2E-se.j1E+1);
//        if ( NULL == ( se.fE1 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
//        if ( NULL == ( se.fE2 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
//        if ( NULL == ( se.rEx = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
//        if ( NULL == ( se.rEz = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
//        for ( size_t i=se.i1E, n=0; i<=se.i2E; ++i ) {
//            for ( size_t j=se.j1E; j<=se.j2E; ++j, ++n ) {
//                size_t ij = (i+g.ab.np)*g.nz2+j+g.ab.np;
//                se.fE1[n] = (1.-mp.phi[ij])*mp.rho_s[ij]*1.e6;
//                se.fE2[n] = mp.phi[ij]*mp.rho_f[ij]*1.e6;
//                se.rEx[n] = (mp.T1[ij]-1.) / (1./mp.phi[ij] - 1.);
//                se.rEz[n] = (mp.T3[ij]-1.) / (1./mp.phi[ij] - 1.);
//            }
//        }
//    }

    //
    // For effective computation:
    // storing   1/rho                                  into   rho
    //           (1-1/eta)                              into   eta
	//           1/tau_s                                into   tau_s
    //
    for ( size_t n=0; n<nnodes; ++n ) {        
        c.rho_i[n] = 1./c.rho_i[n];
        c.rho_j[n] = 1./c.rho_j[n];
        
        c.eta2[n]       = 1.-1./c.eta2[n];
        c.eta2_ij[n]    = 1.-1./c.eta2_ij[n];
		c.tau_s_2[n]    = 1./c.tau_s_2[n];
		c.tau_s_2_ij[n] = 1./c.tau_s_2_ij[n];
    }
	for ( size_t m=0; m<mp.L; ++m ) {
		for ( size_t n=0; n<nnodes; ++n ) {
			c.eta1[m][n] = (1.-1./c.eta1[m][n]) / mp.L;
			c.tau_s_1[m][n] = 1./c.tau_s_1[m][n];
		}
	}

	
	
    free ( mp.K );
    free ( mp.mu );
    free ( mp.rho );
    free ( mp.epsilon );
    free ( mp.delta );
    free ( mp.Q2 );
    free ( mp.f2 );
    for ( size_t m=0; m<mp.L; ++m ) {
        free ( mp.Q1[m] );
        free ( mp.f1[m] );
    }
    free ( mp.Q1 );
    free ( mp.f1 );
    
    for ( size_t n=0; n<(7+mp.L)*nnodes; ++n )
        W[n] = Ws[n] = delta[n] = 0.0;
    
    if ( verbose ) {
        fprintf(stdout, "done.\n");
        fflush(stdout);
    }

	write_field_nc(c.c11, "c11", "MPa", &g, &out, 1);
	write_field_nc(c.c13, "c13", "MPa", &g, &out, 1);
	write_field_nc(c.c33, "c33", "MPa", &g, &out, 1);
	write_field_nc(c.c55, "c55", "MPa", &g, &out, 1);
	write_field_nc(c.K0,  "K0", "MPa", &g, &out, 1);

    
    tau_xx = &(W[0]);
    tau_zz = &(W[nnodes]);
    tau_xz = &(W[2*nnodes]);
    v_x    = &(W[3*nnodes]);
    v_z    = &(W[4*nnodes]);

    init_fftw_data_ve(Wtmp, &fd, &g);
	
	//
	//  Main loop
	//
    
    nsteps = params.duration / c.dt;
    if ( nsteps*c.dt < params.duration ) nsteps++;
    
    if ( verbose )
        fprintf(stdout, "done.\nStarting main loop (%zd iterations)\n", nsteps);

    for ( size_t it=itstart; it<nsteps; ++it ) {
        
		double t = c.dt*(it+1);
		if ( (it+1)%50 == 0 && verbose ) {
            fprintf(stdout, "  Iteration %zd, t = %g s\n", it+1, t); fflush(stdout); }
		
        // update  components of Ws
        for ( size_t n=0; n<(7+mp.L)*nnodes; ++n )
            Ws[n] = W[n];
        
		// update PML memory variables
		for ( size_t n=0; n<2*g.ab.np*g.nz2; ++n ) {
			mem[0].dx_txx[n] = mem[1].dx_txx[n];
			mem[0].dx_txz[n] = mem[1].dx_txz[n];
			mem[0].dx_vx[n]  = mem[1].dx_vx[n];
			mem[0].dx_vz[n]  = mem[1].dx_vz[n];
		}
		for ( size_t n=0; n<2*g.ab.np*g.nx2; ++n ) {
			mem[0].dz_txz[n] = mem[1].dz_txz[n];
			mem[0].dz_tzz[n] = mem[1].dz_tzz[n];
			mem[0].dz_vz[n]  = mem[1].dz_vz[n];
			mem[0].dz_vx[n]  = mem[1].dz_vx[n];
		}
		
        // eq (80) of Carcione (1996)
        for ( size_t n=0; n<(7+mp.L)*nnodes; ++n )
            Wtmp[n] = Ws[n];  // -> using Wtmp rather that Ws for first ride because fftw plans use Wtmp
        
        propagateVE_CPML(Wtmp, delta, &g, &c, &fd, mem, &fp1, mp.L);
        
        // eq (81) of Carcione (1996)
        for ( size_t n=0; n<(7+mp.L)*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt/2. * delta[n];
			W[n]    = Ws[n] + c.dt/6. * delta[n];
		}
		// add source term
		add_src_ve_rk4(&src, W, Wtmp, c.dt/6., c.dt/2., it, nnodes);
		
        propagateVE_CPML(Wtmp, delta, &g, &c, &fd, mem, &fp23, mp.L);
		
        // eq (82) of Carcione (1996)
        for ( size_t n=0; n<(7+mp.L)*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt/2. * delta[n];
			W[n]   +=         c.dt/3. * delta[n];
		}
		// add source term
		add_src_ve_rk4(&src, W, Wtmp, c.dt/3., c.dt/2., it, nnodes);
		
        propagateVE_CPML(Wtmp, delta, &g, &c, &fd, mem, &fp23, mp.L);
        
        // eq (83) of Carcione (1996)
        for ( size_t n=0; n<(7+mp.L)*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt    * delta[n];
			W[n]   +=         c.dt/3. * delta[n];
		}
		// add source term
		add_src_ve_rk4(&src, W, Wtmp, c.dt/3., c.dt, it, nnodes);
		
        propagateVE_CPML(Wtmp, delta, &g, &c, &fd, mem, &fp4, mp.L);
        
        // eq (79) of Carcione (1996)
        for ( size_t n=0; n<(7+mp.L)*nnodes; ++n )
            W[n] += c.dt/6. * delta[n];
		// add source term
		add_src_ve_rk4(&src, W, Wtmp, c.dt/6., 0.0, it, nnodes);
		
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
					data = NULL;
                    break;
            }
			if ( data == NULL ) continue;
			
            if ( ( fabs(dt) < 0.99*c.dt ) ||
                ( ( t > out.r[nr].t0 && out.r[nr].dt > 0.0 && ( fabs(rem) < 1.e-8 ) ) ) ) {
				
                if ( out.r[nr].type == TRACE && params.segy==1 )
                    fill_data_segy(data, t, &out, nr);
                else if ( out.r[nr].type == TRACE && params.segy==0 )
                    write_trace(data, t, &out, nr);
                else if ( out.r[nr].type == SNAPSHOT ) {
					if ( verbose > 1 )
                        printf("Writing snapshot of %s, t = %g s\n", component[out.r[nr].comp], t);
                    write_snapshot_nc(data, t, &g, &out, nr, params.plotStrips);
                }
            }
        }
		
//		if ( params.saveEnergy ) {
//			double E = 0.0;
//			for ( size_t i=se.i1E, n=0; i<=se.i2E; ++i ) {
//				for ( size_t j=se.j1E; j<=se.j2E; ++j, ++n ) {
//					size_t ij = (i+g.ab.np)*g.nz2+j+g.ab.np;
//					
//					// eq. 7.171 of Carcione (2007)  book
//					E += se.fE1[n] * ( v_x[ij]*v_x[ij] + v_z[ij]*v_z[ij] );
//					E -= se.rEx[n] * (v_x[ij]-q_x[ij]) * (v_x[ij]-q_x[ij]);
//					E -= se.rEz[n] * (v_z[ij]-q_z[ij]) * (v_z[ij]-q_z[ij]);
//					E += se.fE2[n] * ( q_x[ij]*q_x[ij] + q_z[ij]*q_z[ij] );
//				}
//			}
//			E *= 0.5;
//			fprintf(se.fid, "%le   %le\n", t, E);
//		}
//		
//		if ( params.checkpoint==1 && it>0 && it%params.chkpt_inc==0 ) {
//			if ( verbose > 1 )
//				printf("Saving checkpoint data at iteration %zd\n", it);
//			save_checkpointVTI(it, &g, &params, &src, &out, &se, &fp1, &fp23, &fp4,
//							   mem, &c, coeff_v_1, coeff_v_3,
//							   coeff_q_1, coeff_q_3, W, Ws, Wtmp, delta);
//		}
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
	
//	if ( params.saveEnergy == 1 ) {
//		fclose( se.fid );
//		free ( se.fE1 );
//		free ( se.fE2 );
//		free ( se.rEx );
//		free ( se.rEz );
//	}
	
	for (size_t n=0; n<src.nsrc; ++n) free( src.s[n].fct );
	free ( src.s );
    free_fftw_data_ve(&fd);
	
    free ( delta );
    free ( Wtmp );
    free ( Ws );
    free ( W );
    for ( size_t m=0; m<mp.L; ++m ) {
        free ( c.eta1[m] );
        free ( c.tau_s_1[m] );
        free ( c.epsilon1[m] );
    }
    free ( c.eta1 );
    free ( c.tau_s_1 );
    free ( c.epsilon1 );
    free ( c.rho_j );
    free ( c.rho_i );
    free ( c.eta2_ij );
    free ( c.eta2 );
    free ( c.tau_s_2_ij );
    free ( c.tau_s_2 );
    free ( c.epsilon3 );
    free ( c.epsilon2 );
    free ( c.c55_0_ij );
    free ( c.c55_0 );
    free ( c.K0 );
	free ( c.c55 );
	free ( c.c33 );
	free ( c.c13 );
	free ( c.c11 );
	
	free_cpml_ve(&fp1, &fp23, &fp4, mem);
	
    if ( verbose ) fprintf(stdout, "\nEnd of computation\n");
    
    return 0;
}

