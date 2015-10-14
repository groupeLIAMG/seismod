//
//  main.c
//  check_vti_model
//
//  Created by Bernard Giroux on 2012-10-11.
//
//

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
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "V.h"
#include "io_utils.h"
#include "pml.h"

int verbose=0;


int main(int argc, char * const argv[])
{

	struct inputParams params;
	struct sourceParams src;
    struct outputParams out;
    struct computationVariablesVTI c;
    struct fac_pml fp1, fp23, fp4;
	struct mem_pml mem[3];
    struct materialPropertiesVTI mp;

    //
    // Grid properties
    //
    struct grid g;
	
	const double pi = 4.0*atan(1.0);
    
	double *m1;        // term in Darcy's law
    double *m3;        // term in Darcy's law
    double *rho;       // composite density                           [ kg/m^3 ]
	
	double *coeff_v_1; // coefficients for analytical solution of the stiff system
    double *coeff_v_3;
    double *coeff_q_1;
    double *coeff_q_3;

	
	params.ab = &g.ab;
	set_defaults(&g, &params);
    process_args(argc, argv, &params);
	read_grid_params(params.modelfile, &g);

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

	alloc_cpml(&fp1, &fp23, &fp4, mem, &g);
	
	read_modelVTI(params.modelfile, &g, &mp);

	//
	// Read in source params
	//
	read_source(params.sourcefile, &g, &src);

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
		c.alpha1[n] = 1. - 1./3. * (c.c11[n]+c12+c.c13[n])/mp.K_s[n];   // eq (9)
		c.alpha3[n] = 1. - 1./3. * (2.*c.c13[n] + c.c33[n])/mp.K_s[n];  // eq (10)
		
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
			c.m3[ij]      = 0.5 * ( m3[ij]       + m3[ijj] );
			
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
    
    
    
    double *veloc_qP1_x, *veloc_qP2_x, *veloc_qS_x, *veloc_qP1_z, *veloc_qP2_z, *veloc_qS_z;
    double *atten_qP1_x, *atten_qP2_x, *atten_qS_x, *atten_qP1_z, *atten_qP2_z, *atten_qS_z;
    if ( NULL == ( veloc_qP1_x  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( veloc_qP2_x  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( veloc_qS_x   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( veloc_qP1_z  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( veloc_qP2_z  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( veloc_qS_z   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( atten_qP1_x  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( atten_qP2_x  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( atten_qS_x   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( atten_qP1_z  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( atten_qP2_z  = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( atten_qS_z   = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	
    double omega = 2.*pi*src.s[0].f;
    double complex V[3];
    for ( size_t n=0; n<nnodes; ++n ) {
        
        // we need undrained elastic components (the diff eqns are solved with variable espilon in carcione99b)
        
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
		
		veloc_qP1_x[n] = creal(V[0]);
		veloc_qP2_x[n] = creal(V[1]);
		veloc_qS_x[n]  = creal(V[2]);
		atten_qP1_x[n] = omega*cimag(1./V[0]);
		atten_qP2_x[n] = omega*cimag(1./V[1]);
		atten_qS_x[n]  = omega*cimag(1./V[2]);
	
        V[0] = VqP1(0., 1., cu11, cu13, cu33, c.c55[n],
                    c.alpha1[n], c.alpha3[n], c.M[n], rho[n], mp.rho_f[n],
                    c.nk1[n], c.nk3[n], m1[n], m3[n], omega);
        V[1] = VqP2(0., 1., cu11, cu13, cu33, c.c55[n],
                    c.alpha1[n], c.alpha3[n], c.M[n], rho[n], mp.rho_f[n],
                    c.nk1[n], c.nk3[n], m1[n], m3[n], omega);
        V[2] = VqS(0., 1., cu11, cu13, cu33, c.c55[n],
                   c.alpha1[n], c.alpha3[n], c.M[n], rho[n], mp.rho_f[n],
                   c.nk1[n], c.nk3[n], m1[n], m3[n], omega);
		
		veloc_qP1_z[n] = creal(V[0]);
		veloc_qP2_z[n] = creal(V[1]);
		veloc_qS_z[n]  = creal(V[2]);
		atten_qP1_z[n] = omega*cimag(1./V[0]);
		atten_qP2_z[n] = omega*cimag(1./V[1]);
		atten_qS_z[n]  = omega*cimag(1./V[2]);
		
    }
    
	strcpy(out.basename, params.basename);

    write_field_nc(veloc_qP1_x, "veloc_qP1_x", "m/s", &g, &out, 1);
    write_field_nc(veloc_qP2_x, "veloc_qP2_x", "m/s", &g, &out, 1);
    write_field_nc(veloc_qS_x, "veloc_qS_x", "m/s", &g, &out, 1);
    write_field_nc(veloc_qP1_z, "veloc_qP1_z", "m/s", &g, &out, 1);
    write_field_nc(veloc_qP2_z, "veloc_qP2_z", "m/s", &g, &out, 1);
    write_field_nc(veloc_qS_z, "veloc_qS_z", "m/s", &g, &out, 1);
    
    write_field_nc(atten_qP1_x, "atten_qP1_x", "1/m", &g, &out, 1);
    write_field_nc(atten_qP2_x, "atten_qP2_x", "1/m", &g, &out, 1);
    write_field_nc(atten_qS_x, "atten_qS_x", "1/m", &g, &out, 1);
    write_field_nc(atten_qP1_z, "atten_qP1_z", "1/m", &g, &out, 1);
    write_field_nc(atten_qP2_z, "atten_qP2_z", "1/m", &g, &out, 1);
    write_field_nc(atten_qS_z, "atten_qS_z", "1/m", &g, &out, 1);
    
	free( veloc_qP1_x );
	free( veloc_qP2_x );
	free( veloc_qS_x );
	free( veloc_qP1_z );
	free( veloc_qP2_z );
	free( veloc_qS_z );
	free( atten_qP1_x );
	free( atten_qP2_x );
	free( atten_qS_x );
	free( atten_qP1_z );
	free( atten_qP2_z );
	free( atten_qS_z );
    
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

    return 0;
}

