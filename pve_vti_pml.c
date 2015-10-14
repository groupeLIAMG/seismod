/*
 *  pve_vti_pml.c
 *
 *  Created by Bernard Giroux on 11-03-08.
 *
 *  PoroViscoElastic wave propagation in 2D VTI anisotropic media, for one 
 *  attenuation mechanism, on a staggered grid
 *
 *
 *  Created by Bernard Giroux on 11-01-06.
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
 *
 */

/*
 
 ------------------------------------------------------------------------------
 
 Usage
 
 pve_vti [options] -p parameter_file.dat
 
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
 0              # segy, save traces in segy format (0 or 1, 0 by default)
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

#include <sys/stat.h>
#include <errno.h>
#include <float.h>
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
#include "V.h"

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
anisotropic_vti      <-  keyword for type of model (isotropic, anisotropic_vti)\n\
K_m(1,1) K_s(1,1) K_f(1,1) phi(1,1) mu(1,1) rho_s(1,1) rho_f(1,1) T3(1,1) eta(1,1) kappa3(1,1) Q(1,1) f0(1,1) epsilon(1,1) delta(1,1) k3k1(1,1) T3T1(1,1)\n\
K_m(1,2) K_s(1,2) K_f(1,2) phi(1,2) mu(1,2) rho_s(1,2) rho_f(1,2) T3(1,2) eta(1,2) kappa3(1,2) Q(1,2) f0(1,2) epsilon(1,2) delta(1,2) k3k1(1,2) T3T1(1,2)\n\
K_m(2,1) K_s(2,1) K_f(2,1) phi(2,1) mu(2,1) rho_s(2,1) rho_f(2,1) T3(2,1) eta(2,1) kappa3(2,1) Q(2,1) f0(2,1) epsilon(2,1) delta(2,1) k3k1(2,1) T3T1(2,1)\n\
K_m(2,2) K_s(2,2) K_f(2,2) phi(2,2) mu(2,2) rho_s(2,2) rho_f(2,2) T3(2,2) eta(2,2) kappa3(2,2) Q(2,2) f0(2,2) epsilon(2,2) delta(2,2) k3k1(2,2) T3T1(2,2)\n\
...\n\
K_m(nx,nz) K_s(nx,nz) K_f(nx,nz) phi(nx,nz) mu(nx,nz) rho_s(nx,nz) rho_f(nx,nz) T3(nx,nz) eta(nx,nz) kappa3(nx,nz) Q(nx,nz) f0(nx,nz) epsilon(nx,nz) delta(nx,nz) k3k1(nx,nz) T3T1(nx,nz)\n\
--- end of file ---\n\
\n\
Variable  Description                                                     Units\n\
\n\
K_m       bulk modulus of drained matrix                                [ GPa ]\n\
K_s       bulk modulus of the solid                                     [ GPa ]\n\
K_f       bulk modulus of the fluid                                     [ GPa ]\n\
phi       porosity                                                        [ - ]\n\
mu        shear modulus of the matrix                                   [ GPa ]\n\
rho_s     solid density                                              [ kg/m^3 ]\n\
rho_f     fluid density                                              [ kg/m^3 ]\n\
T3        tortuosity (along z)                                            [ - ]\n\
eta       fluid viscosity                                                [ cP ]\n\
kappa3    permeability  (along z)                                        [ mD ]\n\
Q         seismic quality factor                                          [ - ]\n\
f0        relaxation frequency                                           [ Hz ]\n\
epsilon   Thomsen anisotropy parameter                                    [ - ]\n\
delta     Thomsen anisotropy parameter                                    [ - ]\n\
k3k1      permeability anisotropy ratio ( kappa_z / kappa_x )             [ - ]\n\
T3T1      tortuosity anisotropy ratio ( T_z / T_x )                       [ - ]\n\
\n\
note: 1 cP = 1e-3 Pa.s ; 1 mD = 9.869233e-10 m^2\n\
\n\
*** Worthy option ***\n\
\n\
Layered models can be input by giving values for just a vertical profile; they\n\
will be duplicated automatically.  For example:\n\
\n\
--- beginning of file ---\n\
2 2                  <-  nx and nz, number of nodes in x & z\n\
1.0 1.0              <-  dx and dz, grid step size in x & z, in meters\n\
0.0 0.0              <-  origin of the grid, in meters\n\
anisotropic_vti      <-  keyword for type of model (isotropic, anisotropic_vti)\n\
K_m(1,1) K_s(1,1) K_f(1,1) phi(1,1) mu(1,1) rho_s(1,1) rho_f(1,1) T3(1,1) eta(1,1) kappa3(1,1) Q(1,1) f0(1,1) epsilon(1,1) delta(1,1) k3k1(1,1) T3T1(1,1)\n\
K_m(1,2) K_s(1,2) K_f(1,2) phi(1,2) mu(1,2) rho_s(1,2) rho_f(1,2) T3(1,2) eta(1,2) kappa3(1,2) Q(1,2) f0(1,2) epsilon(1,2) delta(1,2) k3k1(1,2) T3T1(1,2)\n\
...\n\
K_m(1,nz) K_s(1,nz) K_f(1,nz) phi(1,nz) mu(1,nz) rho_s(1,nz) rho_f(1,nz) T3(1,nz) eta(1,nz) kappa3(1,nz) Q(1,nz) f0(1,nz) epsilon(1,nz) delta(1,nz) k3k1(1,nz) T3T1(1,nz)\n\
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
Components can be Sx, Sz, Sxz, Sf or Bulk\n\
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
Qz             <- particle velocity component of 2nd record\n\
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





int main (int argc, char *argv[]) {
	
    struct inputParams params;
    struct computationVariablesVTI c;
    struct fac_pml fp1, fp23, fp4;
	struct mem_pml mem[3];
    struct materialPropertiesVTI mp;
    struct fftw_data fd;
	struct sourceParams src;
    struct outputParams out;
	struct saveEnergyVTI se;
	
	// -------------------------------------------------------------------------
    //
    // variables
    //
    // -------------------------------------------------------------------------
	
    const double pi = 4.0*atan(1.0);
//    const double L = 1.0;  // number of attenuation mechanisms
    const char *typeSrc[] = { "Sx", "Sy", "Sz", "Sxy", "Sxz", "Syz", "Bulk", "Sf", "Bulk_s" };
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
    
//    double *tau_xxtmp;
//    double *tau_zztmp;
//    double *tau_xztmp;
//    double *ptmp;
	
    //double *e;         // memory variable (one att. mechanism),
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
	
	size_t itstart = 0;
    
	params.ab = &g.ab;
	set_defaults(&g, &params);
    process_args(argc, argv, &params);
	
	if ( params.checkpoint == 1 ) {
		if ( verbose ) {
			fprintf(stdout, "\n *** pve_vti - PoroViscoElastic wave propagation in 2D VTI anisotropic media ***\n\n");
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
			fprintf(stdout, "\n *** pve_vti - PoroViscoElastic wave propagation in 2D VTI anisotropic media ***\n\n");
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
	
	if ( NULL == ( W      = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( Ws     = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( Wtmp   = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( delta  = (double *) malloc( 9*nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	
	alloc_cpml(&fp1, &fp23, &fp4, mem, &g);
	
	if ( params.checkpoint == 1 ) {
		        
		if ( verbose ) {
			fprintf(stdout, "done.\nReading main variables from checkpoint file ... ");
			fflush(stdout);
		}
		read_checkpoint2VTI(params.checkpointfile, &g, &params, &src, &out, &se,
							&fp1, &fp23, &fp4, mem, &c, coeff_v_1, coeff_v_3,
							coeff_q_1, coeff_q_3, W, Ws, Wtmp, delta);
		
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
		read_modelVTI(params.modelfile, &g, &mp);
		
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
				double lambda_s = -(c.nk3[ij]) * beta22;            // eq (71) of Carcione (1996)
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
        
        if ( params.check_model > 0 ) {
            
            for ( size_t n=0; n<nnodes; ++n ) {
                if ( c.alpha1[n] < 0.0 || c.alpha3[n] < 0.0 ) {
                    fprintf(stdout, "\n***\nWarning: bulk modulus of drained matrix larger than modulus of the solid\n***\n");
                    break;
                }
            }
            if ( verbose ) {
                fprintf(stdout, "done.\nComputing phase velocity and attenuation ... ");
                fflush(stdout);
            }
            
            double v_qP1_x_min, v_qP2_x_min, v_qS_x_min, v_qP1_z_min, v_qP2_z_min, v_qS_z_min;
            double v_qP1_x_max, v_qP2_x_max, v_qS_x_max, v_qP1_z_max, v_qP2_z_max, v_qS_z_max;
            double a_qP1_x_min, a_qP2_x_min, a_qS_x_min, a_qP1_z_min, a_qP2_z_min, a_qS_z_min;
            double a_qP1_x_max, a_qP2_x_max, a_qS_x_max, a_qP1_z_max, a_qP2_z_max, a_qS_z_max;
            
            v_qP1_x_min=v_qP2_x_min=v_qS_x_min=v_qP1_z_min=v_qP2_z_min=v_qS_z_min=DBL_MAX;
            v_qP1_x_max=v_qP2_x_max=v_qS_x_max=v_qP1_z_max=v_qP2_z_max=v_qS_z_max=-DBL_MAX;
            a_qP1_x_min=a_qP2_x_min=a_qS_x_min=a_qP1_z_min=a_qP2_z_min=a_qS_z_min=DBL_MAX;
            a_qP1_x_max=a_qP2_x_max=a_qS_x_max=a_qP1_z_max=a_qP2_z_max=a_qS_z_max=-DBL_MAX;
            
            double *vel_qP1_x, *vel_qP2_x, *vel_qS_x, *vel_qP1_z, *vel_qP2_z, *vel_qS_z;
            double *att_qP1_x, *att_qP2_x, *att_qS_x, *att_qP1_z, *att_qP2_z, *att_qS_z;
            if ( NULL == ( vel_qP1_x  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            if ( NULL == ( vel_qP2_x  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            if ( NULL == ( vel_qS_x   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            if ( NULL == ( vel_qP1_z  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            if ( NULL == ( vel_qP2_z  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            if ( NULL == ( vel_qS_z   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            if ( NULL == ( att_qP1_x  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            if ( NULL == ( att_qP2_x  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            if ( NULL == ( att_qS_x   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            if ( NULL == ( att_qP1_z  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            if ( NULL == ( att_qP2_z  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            if ( NULL == ( att_qS_z   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            
            double omega = 2.*pi*f0;
            double complex V[3];
            for ( size_t n=0; n<nnodes; ++n ) {
                
                // we need undrained elastic components (the diff eqns are solved with variable epsilon in carcione99b)
                
                double cu11 = c.c11[n] + c.alpha1[n]*c.alpha1[n]*c.M[n];
                double cu33 = c.c33[n] + c.alpha3[n]*c.alpha3[n]*c.M[n];
                double cu13 = c.c13[n] + c.alpha1[n]*c.alpha3[n]*c.M[n];
                
                V[0] = VqP1(1., 0., cu11, cu13, cu33, c.c55[n],
                            c.alpha1[n], c.alpha3[n], c.M[n], rho[n], mp.rho_f[n],
                            c.nk1[n], c.nk3[n], m1[n], m3[n], omega);
                V[1] = VqP2(1., 0., cu11, cu13, cu33, c.c55[n],
                            c.alpha1[n], c.alpha3[n], c.M[n], rho[n], mp.rho_f[n],
                            c.nk1[n], c.nk3[n], m1[n], m3[n], omega);
                V[2] = VqS(1., 0., cu11, cu13, cu33, c.c55[n],
                           c.alpha1[n], c.alpha3[n], c.M[n], rho[n], mp.rho_f[n],
                           c.nk1[n], c.nk3[n], m1[n], m3[n], omega);
                
                vel_qP1_x[n] = 1./creal(1./V[0]);
                vel_qP2_x[n] = 1./creal(1./V[1]);
                vel_qS_x[n]  = 1./creal(1./V[2]);
                att_qP1_x[n] = -omega*cimag(1./V[0]);
                att_qP2_x[n] = -omega*cimag(1./V[1]);
                att_qS_x[n]  = -omega*cimag(1./V[2]);
                
                v_qP1_x_min = (v_qP1_x_min<vel_qP1_x[n]) ? v_qP1_x_min : vel_qP1_x[n];
                v_qP2_x_min = (v_qP2_x_min<vel_qP2_x[n]) ? v_qP2_x_min : vel_qP2_x[n];
                v_qS_x_min  = ( v_qS_x_min< vel_qS_x[n]) ?  v_qS_x_min :  vel_qS_x[n];

                a_qP1_x_min = (a_qP1_x_min<att_qP1_x[n]) ? a_qP1_x_min : att_qP1_x[n];
                a_qP2_x_min = (a_qP2_x_min<att_qP2_x[n]) ? a_qP2_x_min : att_qP2_x[n];
                a_qS_x_min  = ( a_qS_x_min< att_qS_x[n]) ?  a_qS_x_min :  att_qS_x[n];
                
                v_qP1_x_max = (v_qP1_x_max>vel_qP1_x[n]) ? v_qP1_x_max : vel_qP1_x[n];
                v_qP2_x_max = (v_qP2_x_max>vel_qP2_x[n]) ? v_qP2_x_max : vel_qP2_x[n];
                v_qS_x_max  = ( v_qS_x_max> vel_qS_x[n]) ?  v_qS_x_max :  vel_qS_x[n];
                
                a_qP1_x_max = (a_qP1_x_max>att_qP1_x[n]) ? a_qP1_x_max : att_qP1_x[n];
                a_qP2_x_max = (a_qP2_x_max>att_qP2_x[n]) ? a_qP2_x_max : att_qP2_x[n];
                a_qS_x_max  = ( a_qS_x_max> att_qS_x[n]) ?  a_qS_x_max :  att_qS_x[n];
                
                V[0] = VqP1(0., 1., cu11, cu13, cu33, c.c55[n],
                            c.alpha1[n], c.alpha3[n], c.M[n], rho[n], mp.rho_f[n],
                            c.nk1[n], c.nk3[n], m1[n], m3[n], omega);
                V[1] = VqP2(0., 1., cu11, cu13, cu33, c.c55[n],
                            c.alpha1[n], c.alpha3[n], c.M[n], rho[n], mp.rho_f[n],
                            c.nk1[n], c.nk3[n], m1[n], m3[n], omega);
                V[2] = VqS(0., 1., cu11, cu13, cu33, c.c55[n],
                           c.alpha1[n], c.alpha3[n], c.M[n], rho[n], mp.rho_f[n],
                           c.nk1[n], c.nk3[n], m1[n], m3[n], omega);
                
                vel_qP1_z[n] = 1./creal(1./V[0]);
                vel_qP2_z[n] = 1./creal(1./V[1]);
                vel_qS_z[n]  = 1./creal(1./V[2]);
                att_qP1_z[n] = -omega*cimag(1./V[0]);
                att_qP2_z[n] = -omega*cimag(1./V[1]);
                att_qS_z[n]  = -omega*cimag(1./V[2]);
                
                v_qP1_z_min = (v_qP1_z_min<vel_qP1_z[n]) ? v_qP1_z_min : vel_qP1_z[n];
                v_qP2_z_min = (v_qP2_z_min<vel_qP2_z[n]) ? v_qP2_z_min : vel_qP2_z[n];
                v_qS_z_min  = ( v_qS_z_min< vel_qS_z[n]) ?  v_qS_z_min :  vel_qS_z[n];
                
                a_qP1_z_min = (a_qP1_z_min<att_qP1_z[n]) ? a_qP1_z_min : att_qP1_z[n];
                a_qP2_z_min = (a_qP2_z_min<att_qP2_z[n]) ? a_qP2_z_min : att_qP2_z[n];
                a_qS_z_min  = ( a_qS_z_min< att_qS_z[n]) ?  a_qS_z_min :  att_qS_z[n];
                
                v_qP1_z_max = (v_qP1_z_max>vel_qP1_z[n]) ? v_qP1_z_max : vel_qP1_z[n];
                v_qP2_z_max = (v_qP2_z_max>vel_qP2_z[n]) ? v_qP2_z_max : vel_qP2_z[n];
                v_qS_z_max  = ( v_qS_z_max> vel_qS_z[n]) ?  v_qS_z_max :  vel_qS_z[n];
                
                a_qP1_z_max = (a_qP1_z_max>att_qP1_z[n]) ? a_qP1_z_max : att_qP1_z[n];
                a_qP2_z_max = (a_qP2_z_max>att_qP2_z[n]) ? a_qP2_z_max : att_qP2_z[n];
                a_qS_z_max  = ( a_qS_z_max> att_qS_z[n]) ?  a_qS_z_max :  att_qS_z[n];
                
            }
            if ( verbose ) {
                fprintf(stdout, "done.\n");
                
				fprintf(stdout, "Max/min phase velocities in X:\n");
				fprintf(stdout, "   Src frequency (%lg Hz)   - slow P : %lg    %lg\n", f0, v_qP2_x_max, v_qP2_x_min);
				fprintf(stdout, "                             fast P : %lg    %lg\n", v_qP1_x_max, v_qP1_x_min);
				fprintf(stdout, "                                  S : %lg    %lg\n", v_qS_x_max, v_qS_x_min);
				fprintf(stdout, "Max/min attenuation coefficient in X:\n");
				fprintf(stdout, "   Src frequency (%lg Hz)   - slow P : %lg    %lg\n", f0, a_qP2_x_max, a_qP2_x_min);
				fprintf(stdout, "                             fast P : %lg    %lg\n", a_qP1_x_max, a_qP1_x_min);
				fprintf(stdout, "                                  S : %lg    %lg\n", a_qS_x_max, a_qS_x_min);
				fprintf(stdout, "Max/min phase velocities in Z:\n");
				fprintf(stdout, "   Src frequency (%lg Hz)   - slow P : %lg    %lg\n", f0, v_qP2_z_max, v_qP2_z_min);
				fprintf(stdout, "                             fast P : %lg    %lg\n", v_qP1_z_max, v_qP1_z_min);
				fprintf(stdout, "                                  S : %lg    %lg\n", v_qS_z_max, v_qS_z_min);
				fprintf(stdout, "Max/min attenuation coefficient in Z:\n");
				fprintf(stdout, "   Src frequency (%lg Hz)   - slow P : %lg    %lg\n", f0, a_qP2_z_max, a_qP2_z_min);
				fprintf(stdout, "                             fast P : %lg    %lg\n", a_qP1_z_max, a_qP1_z_min);
				fprintf(stdout, "                                  S : %lg    %lg\n", a_qS_z_max, a_qS_z_min);                
                fflush(stdout);
            }
            
            if ( params.check_model > 1 ) {
                if ( verbose > 1 ) {
                    fprintf(stdout, "  Saving phase velocity and attenuation ...");
                    fflush(stdout);
                }
            
                strcpy(out.basename, params.basename);
                
                write_field_nc(vel_qP1_x, "vel_qP1_x", "m/s", &g, &out, 1);
                write_field_nc(vel_qP2_x, "vel_qP2_x", "m/s", &g, &out, 1);
                write_field_nc(vel_qS_x, "vel_qS_x", "m/s", &g, &out, 1);
                write_field_nc(vel_qP1_z, "vel_qP1_z", "m/s", &g, &out, 1);
                write_field_nc(vel_qP2_z, "vel_qP2_z", "m/s", &g, &out, 1);
                write_field_nc(vel_qS_z, "vel_qS_z", "m/s", &g, &out, 1);
                
                write_field_nc(att_qP1_x, "att_qP1_x", "1/m", &g, &out, 1);
                write_field_nc(att_qP2_x, "att_qP2_x", "1/m", &g, &out, 1);
                write_field_nc(att_qS_x, "att_qS_x", "1/m", &g, &out, 1);
                write_field_nc(att_qP1_z, "att_qP1_z", "1/m", &g, &out, 1);
                write_field_nc(att_qP2_z, "att_qP2_z", "1/m", &g, &out, 1);
                write_field_nc(att_qS_z, "att_qS_z", "1/m", &g, &out, 1);
            }
            
            free( vel_qP1_x );
            free( vel_qP2_x );
            free( vel_qS_x );
            free( vel_qP1_z );
            free( vel_qP2_z );
            free( vel_qS_z );
            free( att_qP1_x );
            free( att_qP2_x );
            free( att_qS_x );
            free( att_qP1_z );
            free( att_qP2_z );
            free( att_qS_z );

        }
		
		if ( params.saveEnergy == 1 ) {
			size_t npts = (se.i2E-se.i1E+1)*(se.j2E-se.j1E+1);
			if ( NULL == ( se.fE1 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
			if ( NULL == ( se.fE2 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
			if ( NULL == ( se.rEx = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
			if ( NULL == ( se.rEz = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
			for ( size_t i=se.i1E, n=0; i<=se.i2E; ++i ) {
				for ( size_t j=se.j1E; j<=se.j2E; ++j, ++n ) {
					size_t ij = (i+g.ab.np)*g.nz2+j+g.ab.np;
					se.fE1[n] = (1.-mp.phi[ij])*mp.rho_s[ij]*1.e6;
					se.fE2[n] = mp.phi[ij]*mp.rho_f[ij]*1.e6;
					se.rEx[n] = (mp.T1[ij]-1.) / (1./mp.phi[ij] - 1.);
					se.rEz[n] = (mp.T3[ij]-1.) / (1./mp.phi[ij] - 1.);
				}
			}
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
		
	    for ( size_t n=0; n<9*nnodes; ++n )
			W[n] = Ws[n] = delta[n] = 0.0;
		
		if ( verbose ) {
			fprintf(stdout, "done.\n");
			fflush(stdout);
		}
		
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
	
//    tau_xxtmp = &(Wtmp[0]);
//    tau_zztmp = &(Wtmp[nnodes]);
//    tau_xztmp = &(Wtmp[2*nnodes]);
//    ptmp      = &(Wtmp[3*nnodes]);
	
    init_fftw_data(Wtmp, &fd, &g);
	
	//
	//  Main loop
	//
    
    nsteps = params.duration / c.dt;
    if ( nsteps*c.dt < params.duration ) nsteps++;
    
    if ( verbose )
        fprintf(stdout, "done.\nStarting main loop (%zd iterations)\n", nsteps);
    for ( size_t it=itstart; it<nsteps; ++it ) {
        
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
        
		// update PML memory variables
		for ( size_t n=0; n<2*g.ab.np*g.nz2; ++n ) {
			mem[0].dx_txx[n] = mem[1].dx_txx[n];
			mem[0].dx_p[n]   = mem[1].dx_p[n];
			mem[0].dx_txz[n] = mem[1].dx_txz[n];
			mem[0].dx_vx[n]  = mem[1].dx_vx[n];
			mem[0].dx_qx[n]  = mem[1].dx_qx[n];
			mem[0].dx_vz[n]  = mem[1].dx_vz[n];
		}
		for ( size_t n=0; n<2*g.ab.np*g.nx2; ++n ) {
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
        
        propagateVTI_CPML(Wtmp, delta, &g, &c, &fd, mem, &fp1);
        
        // eq (81) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt/2. * delta[n];
			W[n]    = Ws[n] + c.dt/6. * delta[n];
		}
		// add source term
		add_src_rk4(&src, W, Wtmp, c.dt/6., c.dt/2., it, nnodes);
		
        propagateVTI_CPML(Wtmp, delta, &g, &c, &fd, mem, &fp23);        
		
        // eq (82) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt/2. * delta[n];
			W[n]   +=         c.dt/3. * delta[n];
		}
		// add source term
		add_src_rk4(&src, W, Wtmp, c.dt/3., c.dt/2., it, nnodes);
		
        propagateVTI_CPML(Wtmp, delta, &g, &c, &fd, mem, &fp23);
        
        // eq (83) of Carcione (1996)
        for ( size_t n=0; n<9*nnodes; ++n ) {
            Wtmp[n] = Ws[n] + c.dt    * delta[n];
			W[n]   +=         c.dt/3. * delta[n];
		}
		// add source term
		add_src_rk4(&src, W, Wtmp, c.dt/3., c.dt, it, nnodes);
		
        propagateVTI_CPML(Wtmp, delta, &g, &c, &fd, mem, &fp4);
        
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
				
                if ( out.r[nr].type == TRACE && params.segy==1 )
                    fill_data_segy(data, t, &out, nr);
                else if ( out.r[nr].type == TRACE && params.segy==0 )
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
					
					// eq. 7.171 of Carcione (2007)  book (Kinetic energy)
					E += se.fE1[n] * ( v_x[ij]*v_x[ij] + v_z[ij]*v_z[ij] );
					E -= se.rEx[n] * (v_x[ij]-q_x[ij]) * (v_x[ij]-q_x[ij]);
					E -= se.rEz[n] * (v_z[ij]-q_z[ij]) * (v_z[ij]-q_z[ij]);
					E += se.fE2[n] * ( q_x[ij]*q_x[ij] + q_z[ij]*q_z[ij] );
				}
			}
			E *= 0.5;
			fprintf(se.fid, "%le   %le\n", t, E);
		}
		
		if ( params.checkpoint==1 && it>0 && it%params.chkpt_inc==0 ) {
			if ( verbose > 1 )
				printf("Saving checkpoint data at iteration %zd\n", it);
			save_checkpointVTI(it, &g, &params, &src, &out, &se, &fp1, &fp23, &fp4,
							   mem, &c, coeff_v_1, coeff_v_3,
							   coeff_q_1, coeff_q_3, W, Ws, Wtmp, delta);
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
	
	if ( params.saveEnergy == 1 ) {
		fclose( se.fid );
		free ( se.fE1 );
		free ( se.fE2 );
		free ( se.rEx );
		free ( se.rEz );
	}
	
	for (size_t n=0; n<src.nsrc; ++n) free( src.s[n].fct );
	free ( src.s );
    free_fftw_data(&fd);
	
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
	
	free_cpml(&fp1, &fp23, &fp4, mem);
	
    if ( verbose ) fprintf(stdout, "\nEnd of computation\n");
    
	return 0;
}