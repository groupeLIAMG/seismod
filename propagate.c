/*
 *  propagate.c
 *
 *  Created by Bernard Giroux on 10-08-29.
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


#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "propagate.h"

extern int verbose;

void init_fftw_data(double *W, struct fftw_data *d, const struct grid *g) {
    
    size_t nnodes = g->nx2 * g->nz2;
    size_t nx = g->nx2;
    size_t nz = g->nz2;
    
    double *tau_xx = &(W[0]);
    double *tau_zz = &(W[nnodes]);
    double *tau_xz = &(W[2*nnodes]);
    double *p      = &(W[3*nnodes]);
    double *v_x    = &(W[4*nnodes]);
    double *v_z    = &(W[5*nnodes]);
    double *q_x    = &(W[6*nnodes]);
    double *q_z    = &(W[7*nnodes]);
    
#ifdef _OPENMP
    if ( verbose >= 1 ) {
        fprintf(stdout, "Initializing fftw data - using %d threads (this may take some time) ... ", omp_get_max_threads());
        fflush(stdout);
    }
    fftw_init_threads();
    fftw_plan_with_nthreads( omp_get_max_threads() );
#else
    if ( verbose >= 1 ) {
        fprintf(stdout, "Initializing fftw data (this may take some time) ... ");
        fflush(stdout);
    }	
#endif

    if ( NULL == ( d->o_x  = (fftw_complex*)    fftw_malloc((nx/2+1)*nz*sizeof(fftw_complex)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->o_z  = (fftw_complex*)    fftw_malloc(nx*(nz/2+1)*sizeof(fftw_complex)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->t1   = (double *)         fftw_malloc(nnodes * sizeof(double))))          { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->t2   = (double *)         fftw_malloc(nnodes * sizeof(double))))          { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->t3   = (double *)         fftw_malloc(nnodes * sizeof(double))))          { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kx_f = (complex double *) fftw_malloc( (nx/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kx_b = (complex double *) fftw_malloc( (nx/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kz_f = (complex double *) fftw_malloc( (nz/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kz_b = (complex double *) fftw_malloc( (nz/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }

	double dn = 8.0*atan(1.0) / (g->dx*nx*nx);
    for (size_t n=0; n<nx/2+1; ++n) {
        double k = n*dn;
        d->kx_f[n] = I*k*cexp(  0.5*I*k*g->dx*nx );
        d->kx_b[n] = I*k*cexp( -0.5*I*k*g->dx*nx );
	}
	dn = 8.0*atan(1.0) / (g->dz*nz*nz);
    for (size_t n=0; n<nz/2+1; ++n) {
        double k = n*dn;
        d->kz_f[n] = I*k*cexp(  0.5*I*k*g->dz*nz );
        d->kz_b[n] = I*k*cexp( -0.5*I*k*g->dz*nz );
	}
    
    char hostname[100], wisdomfile[200];
    gethostname( hostname, 100 );
    sprintf(wisdomfile, "%s/.fftw/%s.wisdom", getenv("HOME"), hostname);
    FILE *wfile = fopen(wisdomfile, "r");
    if ( wfile != NULL ) {
        fftw_import_wisdom_from_file( wfile );
		fclose( wfile );
	}
    
    int n[] = { (int)nx };
    int howmany = (int)nz;
    int inembed[] = { (int)nnodes };
    int istride = (int)nz;
    int idist = 1;
    int onembed[] = { (int)((nx/2+1)*nz) };
    int ostride = 1;
    int odist = nx/2+1;
    
	if ( NULL == ( 
    d->tauxx_x_f = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          tau_xx, inembed, istride, idist,
                                          d->o_x, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->tauxz_x_f = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          tau_xz, inembed, istride, idist,
                                          d->o_x, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->p_x_f     = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          p,      inembed, istride, idist,
                                          d->o_x, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->vx_x_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          v_x,    inembed, istride, idist,
                                          d->o_x, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->vz_x_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          v_z,    inembed, istride, idist,
                                          d->o_x, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->qx_x_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          q_x,    inembed, istride, idist,
                                          d->o_x, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }


	if ( NULL == ( 
    d->tauxx_x_i = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_x, onembed, ostride, odist,
                                          d->t1,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->tauxz_x_i = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_x, onembed, ostride, odist,
                                          d->t1,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->p_x_i     = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_x, onembed, ostride, odist,
                                          d->t2,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }

	if ( NULL == ( 
	d->vx_x_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_x, onembed, ostride, odist,
                                          d->t1,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->vz_x_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_x, onembed, ostride, odist,
                                          d->t2,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->qx_x_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_x, onembed, ostride, odist,
                                          d->t3,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	
    
    n[0] = (int)nz;
    howmany = (int)nx;
    inembed[0] = (int)nnodes;
    istride = 1;
    idist = (int)nz;
    onembed[0] = (int)(nx*(nz/2+1));
    ostride = 1;
    odist = nz/2+1;
    
	if ( NULL == ( 
    d->tauxz_z_f = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          tau_xz, inembed, istride, idist,
                                          d->o_z, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->tauzz_z_f = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          tau_zz, inembed, istride, idist,
                                          d->o_z, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->p_z_f     = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          p,      inembed, istride, idist,
                                          d->o_z, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->vx_z_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          v_x,    inembed, istride, idist,
                                          d->o_z, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->vz_z_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          v_z,    inembed, istride, idist,
                                          d->o_z, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->qz_z_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                          q_z,    inembed, istride, idist,
                                          d->o_z, onembed, ostride, odist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }

	if ( NULL == ( 
    d->tauxz_z_i = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_z, onembed, ostride, odist,
                                          d->t2,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->tauzz_z_i = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_z, onembed, ostride, odist,
                                          d->t2,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->p_z_i     = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_z, onembed, ostride, odist,
                                          d->t2,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->vx_z_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_z, onembed, ostride, odist,
                                          d->t1,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->vz_z_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_z, onembed, ostride, odist,
                                          d->t2,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( 
    d->qz_z_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                          d->o_z, onembed, ostride, odist,
                                          d->t3,  inembed, istride, idist,
                                          FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }

    
    if ( (wfile = fopen(wisdomfile, "w")) != NULL ) {
		fftw_export_wisdom_to_file( wfile );
		fclose( wfile );
	}
}

void free_fftw_data(struct fftw_data *d) {
    
    fftw_destroy_plan( d->tauxx_x_f );
    fftw_destroy_plan( d->tauzz_z_f );
    fftw_destroy_plan( d->tauxz_x_f );
    fftw_destroy_plan( d->tauxz_z_f );
    fftw_destroy_plan( d->p_x_f );
    fftw_destroy_plan( d->p_z_f );
    fftw_destroy_plan( d->vx_x_f );
    fftw_destroy_plan( d->vx_z_f );
    fftw_destroy_plan( d->vz_x_f );
    fftw_destroy_plan( d->vz_z_f );
    fftw_destroy_plan( d->qx_x_f );
    fftw_destroy_plan( d->qz_z_f );

    fftw_destroy_plan( d->tauxx_x_i );
    fftw_destroy_plan( d->tauzz_z_i );
    fftw_destroy_plan( d->tauxz_x_i );
    fftw_destroy_plan( d->tauxz_z_i );
    fftw_destroy_plan( d->p_x_i );
    fftw_destroy_plan( d->p_z_i );
    fftw_destroy_plan( d->vx_x_i );
    fftw_destroy_plan( d->vx_z_i );
    fftw_destroy_plan( d->vz_x_i );
    fftw_destroy_plan( d->vz_z_i );
    fftw_destroy_plan( d->qx_x_i );
    fftw_destroy_plan( d->qz_z_i );
    
    fftw_free( d->o_x );
    fftw_free( d->o_z );
    fftw_free( d->t1 );
    fftw_free( d->t2 );
    fftw_free( d->t3 );
    fftw_free( d->kx_f );
    fftw_free( d->kx_b );
    fftw_free( d->kz_f );
    fftw_free( d->kz_b );
}

void init_fftw_data_ve(double *W, struct fftw_data *d, const struct grid *g) {
    
    size_t nnodes = g->nx2 * g->nz2;
    size_t nx = g->nx2;
    size_t nz = g->nz2;
    
    double *tau_xx = &(W[0]);
    double *tau_zz = &(W[nnodes]);
    double *tau_xz = &(W[2*nnodes]);
    double *v_x    = &(W[3*nnodes]);
    double *v_z    = &(W[4*nnodes]);
    
#ifdef _OPENMP
    if ( verbose >= 1 ) {
        fprintf(stdout, "Initializing fftw data - using %d threads (this may take some time) ... ", omp_get_max_threads());
        fflush(stdout);
    }
    fftw_init_threads();
    fftw_plan_with_nthreads( omp_get_max_threads() );
#else
    if ( verbose >= 1 ) {
        fprintf(stdout, "Initializing fftw data (this may take some time) ... ");
        fflush(stdout);
    }
#endif
    
    if ( NULL == ( d->o_x  = (fftw_complex*)    fftw_malloc((nx/2+1)*nz*sizeof(fftw_complex)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->o_z  = (fftw_complex*)    fftw_malloc(nx*(nz/2+1)*sizeof(fftw_complex)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->t1   = (double *)         fftw_malloc(nnodes * sizeof(double))))          { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->t2   = (double *)         fftw_malloc(nnodes * sizeof(double))))          { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kx_f = (complex double *) fftw_malloc( (nx/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kx_b = (complex double *) fftw_malloc( (nx/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kz_f = (complex double *) fftw_malloc( (nz/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kz_b = (complex double *) fftw_malloc( (nz/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
	double dn = 8.0*atan(1.0) / (g->dx*nx*nx);
    for (size_t n=0; n<nx/2+1; ++n) {
        double k = n*dn;
        d->kx_f[n] = I*k*cexp(  0.5*I*k*g->dx*nx );
        d->kx_b[n] = I*k*cexp( -0.5*I*k*g->dx*nx );
	}
	dn = 8.0*atan(1.0) / (g->dz*nz*nz);
    for (size_t n=0; n<nz/2+1; ++n) {
        double k = n*dn;
        d->kz_f[n] = I*k*cexp(  0.5*I*k*g->dz*nz );
        d->kz_b[n] = I*k*cexp( -0.5*I*k*g->dz*nz );
	}
    
    char hostname[100], wisdomfile[200];
    gethostname( hostname, 100 );
    sprintf(wisdomfile, "%s/.fftw/%s.wisdom", getenv("HOME"), hostname);
    FILE *wfile = fopen(wisdomfile, "r");
    if ( wfile != NULL ) {
        fftw_import_wisdom_from_file( wfile );
		fclose( wfile );
	}
    
    int n[] = { (int)nx };
    int howmany = (int)nz;
    int inembed[] = { (int)nnodes };
    int istride = (int)nz;
    int idist = 1;
    int onembed[] = { (int)((nx/2+1)*nz) };
    int ostride = 1;
    int odist = nx/2+1;
    
	if ( NULL == (
                  d->tauxx_x_f = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        tau_xx, inembed, istride, idist,
                                                        d->o_x, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->tauxz_x_f = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        tau_xz, inembed, istride, idist,
                                                        d->o_x, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->vx_x_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        v_x,    inembed, istride, idist,
                                                        d->o_x, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->vz_x_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        v_z,    inembed, istride, idist,
                                                        d->o_x, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
    
	if ( NULL == (
                  d->tauxx_x_i = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_x, onembed, ostride, odist,
                                                        d->t1,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->tauxz_x_i = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_x, onembed, ostride, odist,
                                                        d->t1,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
	if ( NULL == (
                  d->vx_x_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_x, onembed, ostride, odist,
                                                        d->t1,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->vz_x_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_x, onembed, ostride, odist,
                                                        d->t2,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	
    
    n[0] = (int)nz;
    howmany = (int)nx;
    inembed[0] = (int)nnodes;
    istride = 1;
    idist = (int)nz;
    onembed[0] = (int)(nx*(nz/2+1));
    ostride = 1;
    odist = nz/2+1;
    
	if ( NULL == (
                  d->tauxz_z_f = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        tau_xz, inembed, istride, idist,
                                                        d->o_z, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->tauzz_z_f = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        tau_zz, inembed, istride, idist,
                                                        d->o_z, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->vx_z_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        v_x,    inembed, istride, idist,
                                                        d->o_z, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->vz_z_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        v_z,    inembed, istride, idist,
                                                        d->o_z, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
	if ( NULL == (
                  d->tauxz_z_i = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_z, onembed, ostride, odist,
                                                        d->t2,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->tauzz_z_i = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_z, onembed, ostride, odist,
                                                        d->t2,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->vx_z_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_z, onembed, ostride, odist,
                                                        d->t1,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->vz_z_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_z, onembed, ostride, odist,
                                                        d->t2,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
    
    if ( (wfile = fopen(wisdomfile, "w")) != NULL ) {
		fftw_export_wisdom_to_file( wfile );
		fclose( wfile );
	}
}

void free_fftw_data_ve(struct fftw_data *d) {
    
    fftw_destroy_plan( d->tauxx_x_f );
    fftw_destroy_plan( d->tauzz_z_f );
    fftw_destroy_plan( d->tauxz_x_f );
    fftw_destroy_plan( d->tauxz_z_f );
    fftw_destroy_plan( d->vx_x_f );
    fftw_destroy_plan( d->vx_z_f );
    fftw_destroy_plan( d->vz_x_f );
    fftw_destroy_plan( d->vz_z_f );
    
    fftw_destroy_plan( d->tauxx_x_i );
    fftw_destroy_plan( d->tauzz_z_i );
    fftw_destroy_plan( d->tauxz_x_i );
    fftw_destroy_plan( d->tauxz_z_i );
    fftw_destroy_plan( d->vx_x_i );
    fftw_destroy_plan( d->vx_z_i );
    fftw_destroy_plan( d->vz_x_i );
    fftw_destroy_plan( d->vz_z_i );
    
    fftw_free( d->o_x );
    fftw_free( d->o_z );
    fftw_free( d->t1 );
    fftw_free( d->t2 );
    fftw_free( d->kx_f );
    fftw_free( d->kx_b );
    fftw_free( d->kz_f );
    fftw_free( d->kz_b );
}

void init_fftw_data_ve_sh(double *W, struct fftw_data_ve_sh *d, const struct grid *g) {
    
    size_t nnodes = g->nx2 * g->nz2;
    size_t nx = g->nx2;
    size_t nz = g->nz2;
    
    double *tau_xy = &(W[0]);
    double *tau_yz = &(W[nnodes]);
    double *v_y    = &(W[2*nnodes]);
    
#ifdef _OPENMP
    if ( verbose >= 1 ) {
        fprintf(stdout, "Initializing fftw data - using %d threads (this may take some time) ... ", omp_get_max_threads());
        fflush(stdout);
    }
    fftw_init_threads();
    fftw_plan_with_nthreads( omp_get_max_threads() );
#else
    if ( verbose >= 1 ) {
        fprintf(stdout, "Initializing fftw data (this may take some time) ... ");
        fflush(stdout);
    }
#endif
    
    if ( NULL == ( d->o_x  = (fftw_complex*)    fftw_malloc((nx/2+1)*nz*sizeof(fftw_complex)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->o_z  = (fftw_complex*)    fftw_malloc(nx*(nz/2+1)*sizeof(fftw_complex)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->t1   = (double *)         fftw_malloc(nnodes * sizeof(double))))          { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->t2   = (double *)         fftw_malloc(nnodes * sizeof(double))))          { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kx_f = (complex double *) fftw_malloc( (nx/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kx_b = (complex double *) fftw_malloc( (nx/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kz_f = (complex double *) fftw_malloc( (nz/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( d->kz_b = (complex double *) fftw_malloc( (nz/2+1)*sizeof(complex double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
	double dn = 8.0*atan(1.0) / (g->dx*nx*nx);
    for (size_t n=0; n<nx/2+1; ++n) {
        double k = n*dn;
        d->kx_f[n] = I*k*cexp(  0.5*I*k*g->dx*nx );
        d->kx_b[n] = I*k*cexp( -0.5*I*k*g->dx*nx );
	}
	dn = 8.0*atan(1.0) / (g->dz*nz*nz);
    for (size_t n=0; n<nz/2+1; ++n) {
        double k = n*dn;
        d->kz_f[n] = I*k*cexp(  0.5*I*k*g->dz*nz );
        d->kz_b[n] = I*k*cexp( -0.5*I*k*g->dz*nz );
	}
    
    char hostname[100], wisdomfile[200];
    gethostname( hostname, 100 );
    sprintf(wisdomfile, "%s/.fftw/%s.wisdom", getenv("HOME"), hostname);
    FILE *wfile = fopen(wisdomfile, "r");
    if ( wfile != NULL ) {
        fftw_import_wisdom_from_file( wfile );
		fclose( wfile );
	}
    
    int n[] = { (int)nx };
    int howmany = (int)nz;
    int inembed[] = { (int)nnodes };
    int istride = (int)nz;
    int idist = 1;
    int onembed[] = { (int)((nx/2+1)*nz) };
    int ostride = 1;
    int odist = nx/2+1;
    
	if ( NULL == (
                  d->tauxy_x_f = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        tau_xy, inembed, istride, idist,
                                                        d->o_x, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->vy_x_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        v_y,    inembed, istride, idist,
                                                        d->o_x, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
    
	if ( NULL == (
                  d->tauxy_x_i = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_x, onembed, ostride, odist,
                                                        d->t1,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->vy_x_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_x, onembed, ostride, odist,
                                                        d->t1,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	
    
    n[0] = (int)nz;
    howmany = (int)nx;
    inembed[0] = (int)nnodes;
    istride = 1;
    idist = (int)nz;
    onembed[0] = (int)(nx*(nz/2+1));
    ostride = 1;
    odist = nz/2+1;
    
	if ( NULL == (
                  d->tauyz_z_f = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        tau_yz, inembed, istride, idist,
                                                        d->o_z, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->vy_z_f    = fftw_plan_many_dft_r2c( 1, n, howmany,
                                                        v_y,    inembed, istride, idist,
                                                        d->o_z, onembed, ostride, odist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
	if ( NULL == (
                  d->tauyz_z_i = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_z, onembed, ostride, odist,
                                                        d->t2,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == (
                  d->vy_z_i    = fftw_plan_many_dft_c2r( 1, n, howmany,
                                                        d->o_z, onembed, ostride, odist,
                                                        d->t1,  inembed, istride, idist,
                                                        FFTW_PATIENT)
				  )) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
    
    if ( (wfile = fopen(wisdomfile, "w")) != NULL ) {
		fftw_export_wisdom_to_file( wfile );
		fclose( wfile );
	}
}

void free_fftw_data_ve_sh(struct fftw_data_ve_sh *d) {
    
    fftw_destroy_plan( d->tauxy_x_f );
    fftw_destroy_plan( d->tauyz_z_f );
    fftw_destroy_plan( d->vy_x_f );
    fftw_destroy_plan( d->vy_z_f );
    
    fftw_destroy_plan( d->tauxy_x_i );
    fftw_destroy_plan( d->tauyz_z_i );
    fftw_destroy_plan( d->vy_x_i );
    fftw_destroy_plan( d->vy_z_i );
    
    fftw_free( d->o_x );
    fftw_free( d->o_z );
    fftw_free( d->t1 );
    fftw_free( d->t2 );
    fftw_free( d->kx_f );
    fftw_free( d->kx_b );
    fftw_free( d->kz_f );
    fftw_free( d->kz_b );
}


void propagate(double *W, double *delta, const struct grid *g,
               const struct computationVariables *c, struct fftw_data *d) {
    
    size_t nnodes = g->nx2 * g->nz2;

    //double *q_x     = &(W[6*nnodes]);
    //double *q_z     = &(W[7*nnodes]);
    double *e       = &(W[8*nnodes]);
    
    double *tau2_xx = &(delta[0]);
    double *tau2_zz = &(delta[nnodes]);
    double *tau2_xz = &(delta[2*nnodes]);
    double *p2      = &(delta[3*nnodes]);
    double *v2_x    = &(delta[4*nnodes]);
    double *v2_z    = &(delta[5*nnodes]);
    double *q2_x    = &(delta[6*nnodes]);
    double *q2_z    = &(delta[7*nnodes]);
    double *e2      = &(delta[8*nnodes]);
    
    //
    // Solving system (26) in Carcione and Helle (1999)
    //

	// -------------------------------------------------------------------------
    //
    //  Solve for 
    //
    // D_x^+\tau_{xx} + D_z^-\tau_{xz} = \rho\frac{\partial v_x}{\partial t} + \rho_f\frac{\partial q_x}{\partial t}
    //    and
    // -D_x^+p = \rho_f\frac{\partial v_x}{\partial t} + m\frac{\partial q_x}{\partial t} + \frac{\eta}{\kappa}q_x
    //
    //  Note:
    //           rho_f^2/rho - m                      is stored into   m
    //           rho_f/rho                            is stored into   rho_f
    //           1/rho                                is stored into   rho

    
    // D_x+ tau_xx
    fftw_execute( d->tauxx_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->tauxx_x_i );  // d->t1 now holds D_x tau_xx
    // Dz- tau_xz
    fftw_execute( d->tauxz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->tauxz_z_i );  // d->t2 now holds D_z tau_xz
    
    for ( size_t n=0; n<nnodes; ++n )
        d->t1[n] += d->t2[n];  // d->t1 now holds D_x tau_xx  +  D_z tau_xz
    
    // D_x+ p
    fftw_execute( d->p_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->p_x_i );  // d->t2 now holds D_x p
    
    for ( size_t n=0; n<nnodes; ++n ) {
//        q2_x[n] = (d->t2[n] + c->rho_f_i[n]*d->t1[n] + c->nk_i[n]*q_x[n]) / c->m_i[n];
        q2_x[n] = (d->t2[n] + c->rho_f_i[n]*d->t1[n]) / c->m_i[n];
        v2_x[n] = c->rho_i[n]*d->t1[n] - c->rho_f_i[n]*q2_x[n];
    }    
    
    
    
    // -------------------------------------------------------------------------
    //
    //  Solve for 
    //
    // D_x^+\tau_{xz} + D_z^-\tau_{zz} = \rho\frac{\partial v_z}{\partial t} + \rho_f\frac{\partial q_z}{\partial t}
    //    and
    // -D_z^+p = \rho_f\frac{\partial v_z}{\partial t} + m\frac{\partial q_z}{\partial t} + \frac{\eta}{\kappa}q_z
    
    // D_x- tau_xz
    fftw_execute( d->tauxz_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->tauxz_x_i );  // d->t1 now holds D_x tau_xz
    // D_z+ tau_zz
    fftw_execute( d->tauzz_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->tauzz_z_i );  // d->t2 now holds D_z tau_zz
    
    for ( size_t n=0; n<nnodes; ++n )
        d->t1[n] += d->t2[n];  // d->t1 now holds D_x tau_xz  +  D_z tau_zz
	
    // D_z+ p
    fftw_execute( d->p_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->p_z_i );  // d->t2 now holds D_z p
	
    for ( size_t n=0; n<nnodes; ++n ) {
//        q2_z[n] = (d->t2[n] + c->rho_f_j[n]*d->t1[n] + c->nk_j[n]*q_z[n]) / c->m_j[n];
        q2_z[n] = (d->t2[n] + c->rho_f_j[n]*d->t1[n]) / c->m_j[n];
        v2_z[n] = c->rho_j[n]*d->t1[n] - c->rho_f_j[n]*q2_z[n];
    }    
	
    

    // -------------------------------------------------------------------------
    // \epsilon = \alpha\left(D_x^-v_x + D_z^-v_z\right) + D_x^-q_x + D_z^-q_z
    //
    
    // D_x- v_x
    fftw_execute( d->vx_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->vx_x_i );  // d->t1 now holds D_xv_x
    
    // D_z- v_z
    fftw_execute( d->vz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->vz_z_i );  // d->t2 now holds D_zv_z
    
    for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] = c->alpha[n] * ( d->t1[n]+d->t2[n] );
    
    // D_x- q_x
    fftw_execute( d->qx_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->qx_x_i );  // d->t3 now holds D_xq_x
    for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] += d->t3[n];
    
    // D_z- q_z
    fftw_execute( d->qz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->qz_z_i );  // d->t3 now holds D_zq_z
    for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] += d->t3[n];
    
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial e}{\partial t} = -\frac{1}{\tau_\sigma}\left[M\left(L+\varphi\right)^{-1}\varphi\epsilon+e\right]
    //
    //      1/tau_s                              is stored into   tau_s
    //      M*varphi/(1+varphi)                  is stored into   varphi

    for ( size_t n=0; n<nnodes; ++n ) {
        e2[n] = -c->tau_s[n] * (c->varphi[n]*c->epsilon[n] + e[n]);
        //  store   M*epsilon + e  in  epsilon
        c->epsilon[n] = c->M[n]*c->epsilon[n] + e[n];
    }

    // -------------------------------------------------------------------------
    //
    // \frac{\partial p}{\partial t} = -\left(M\epsilon+e\right)+s_f
    //
    for ( size_t n=0; n<nnodes; ++n )
        p2[n] = - c->epsilon[n];
	
	
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{xx}}{\partial t} = E D_x^-v_x + (E-2\mu)D_z^-v_z+\alpha(M\epsilon+e)+s_{x}
    //
    // c->E  holds in fact E
    // c->mu holds in fact (E-2mu)
    for ( size_t n=0; n<nnodes; ++n )
        tau2_xx[n] = c->E[n]*d->t1[n] + c->mu[n]*d->t2[n] + c->alpha[n]*c->epsilon[n];
	
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{zz}}{\partial t} = (E-2\mu) D_x^-v_x + E D_z^-v_z+\alpha(M\epsilon+e)+s_{z}
    //
    for ( size_t n=0; n<nnodes; ++n )
        tau2_zz[n] = c->mu[n]*d->t1[n] + c->E[n]*d->t2[n] + c->alpha[n]*c->epsilon[n];
	
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{xz}}{\partial t} = \mu\left(D_z^+v_x + D_x^+v_z\right)+s_{xz}
    //
    //    with \mu interpolated at (i+1/2,j+1/2)
    
    // D_z+ v_x
    fftw_execute( d->vx_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->vx_z_i );  // d->t1 now holds D_zv_x

    // D_x+ v_z
    fftw_execute( d->vz_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->vz_x_i );  // d->t2 now holds D_xv_z

    for ( size_t n=0; n<nnodes; ++n )
        tau2_xz[n] = c->mu_ij[n] * ( d->t1[n]+d->t2[n] );
    
}

void propagateCPML(double *W, double *delta, const struct grid *g,
                   const struct computationVariables *c, struct fftw_data *d,
                   struct mem_pml *m, struct fac_pml *p) {
    
    size_t nnodes = g->nx2 * g->nz2;
    size_t npml = g->ab.np;
    
//    double *q_x     = &(W[6*nnodes]);
//    double *q_z     = &(W[7*nnodes]);
    double *e       = &(W[8*nnodes]);
    
    double *tau2_xx = &(delta[0]);
    double *tau2_zz = &(delta[nnodes]);
    double *tau2_xz = &(delta[2*nnodes]);
    double *p2      = &(delta[3*nnodes]);
    double *v2_x    = &(delta[4*nnodes]);
    double *v2_z    = &(delta[5*nnodes]);
    double *q2_x    = &(delta[6*nnodes]);
    double *q2_z    = &(delta[7*nnodes]);
    double *e2      = &(delta[8*nnodes]);
    
    //
    // Solving system (26) in Carcione and Helle (1999)
    //
    
	// -------------------------------------------------------------------------
    //
    //  Solve for 
    //
    // D_x^+\tau_{xx} + D_z^-\tau_{xz} = \rho\frac{\partial v_x}{\partial t} + \rho_f\frac{\partial q_x}{\partial t}
    //    and
    // -D_x^+p = \rho_f\frac{\partial v_x}{\partial t} + m\frac{\partial q_x}{\partial t} + \frac{\eta}{\kappa}q_x
    //
    //  Note:
    //           rho_f^2/rho - m                      is stored into   m
    //           rho_f/rho                            is stored into   rho_f
    //           1/rho                                is stored into   rho
    
    
    // D_x+ tau_xx
    fftw_execute( d->tauxx_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->tauxx_x_i );  // d->t1 now holds D_x tau_xx
    // Dz- tau_xz
    fftw_execute( d->tauxz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->tauxz_z_i );  // d->t2 now holds D_z tau_xz
    
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j
	size_t ipml = 0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_txx[ipml] = p->bh_x[nl]*m[0].dx_txx[ipml] + p->ch_x[nl]*d->t1[ind];  // eq 32 of Martin et al. 2010, with sign corrected before c_x
            d->t1[ind] = d->t1[ind]/p->kh_x[nl] + m[2].dx_txx[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // right
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_txx[ipml] = p->bH_x[nl]*m[0].dx_txx[ipml] + p->cH_x[nl] *d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kH_x[nl] + m[2].dx_txx[ipml];
        }
    }
    // derivative evaluated at i+1/2,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_txz[ipml] = p->b_z[nl]*m[0].dz_txz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_txz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_txz[ipml] = p->b_z[nl]*m[0].dz_txz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_txz[ipml];
        }
    }
    
    for ( size_t n=0; n<nnodes; ++n )
        d->t1[n] += d->t2[n];  // d->t1 now holds D_x tau_xx  +  D_z tau_xz
        
    // D_x+ p
    fftw_execute( d->p_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->p_x_i );  // d->t2 now holds D_x p
    
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_p[ipml] = p->bh_x[nl]*m[0].dx_p[ipml] + p->ch_x[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kh_x[nl] + m[2].dx_p[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // right
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_p[ipml] = p->bH_x[nl]*m[0].dx_p[ipml] + p->cH_x[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kH_x[nl] + m[2].dx_p[ipml];            
        }
    }
    
	for (size_t i=1; i<g->nx2-1; ++i) {
		for ( size_t j=1; j<g->nz2-1; ++j) {
			size_t n = i*g->nz2+j;
//			q2_x[n] = (d->t2[n] + c->rho_f_i[n]*d->t1[n] + c->nk_i[n]*q_x[n]) / c->m_i[n];
			q2_x[n] = (d->t2[n] + c->rho_f_i[n]*d->t1[n]) / c->m_i[n];
			v2_x[n] = c->rho_i[n]*d->t1[n] - c->rho_f_i[n]*q2_x[n];
		}
	}
    
    
    
    // -------------------------------------------------------------------------
    //
    //  Solve for 
    //
    // D_x^+\tau_{xz} + D_z^-\tau_{zz} = \rho\frac{\partial v_z}{\partial t} + \rho_f\frac{\partial q_z}{\partial t}
    //    and
    // -D_z^+p = \rho_f\frac{\partial v_z}{\partial t} + m\frac{\partial q_z}{\partial t} + \frac{\eta}{\kappa}q_z
    
    // D_x- tau_xz
    fftw_execute( d->tauxz_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->tauxz_x_i );  // d->t1 now holds D_x tau_xz
    // D_z+ tau_zz
    fftw_execute( d->tauzz_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->tauzz_z_i );  // d->t2 now holds D_z tau_zz
    
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_txz[ipml] = p->b_x[nl]*m[0].dx_txz[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_txz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {  // right
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_txz[ipml] = p->b_x[nl]*m[0].dx_txz[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_txz[ipml];
        }
    }
    // derivative evaluated at i,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_tzz[ipml] = p->bh_z[nl]*m[0].dz_tzz[ipml] + p->ch_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kh_z[nl] + m[2].dz_tzz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_tzz[ipml] = p->bH_z[nl]*m[0].dz_tzz[ipml] + p->cH_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kH_z[nl] + m[2].dz_tzz[ipml];
        }
    }
    
    for ( size_t n=0; n<nnodes; ++n )
        d->t1[n] += d->t2[n];  // d->t1 now holds D_x tau_xz  +  D_z tau_zz
	
    // D_z+ p
    fftw_execute( d->p_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->p_z_i );  // d->t2 now holds D_z p
	
    // derivative evaluated at i,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_p[ipml] = p->bh_z[nl]*m[0].dz_p[ipml] + p->ch_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kh_z[nl] + m[2].dz_p[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_p[ipml] = p->bH_z[nl]*m[0].dz_p[ipml] + p->cH_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kH_z[nl] + m[2].dz_p[ipml];
        }
    }

	for (size_t i=1; i<g->nx2-1; ++i) {
		for ( size_t j=1; j<g->nz2-1; ++j) {
			size_t n = i*g->nz2+j;
//			q2_z[n] = (d->t2[n] + c->rho_f_j[n]*d->t1[n] + c->nk_j[n]*q_z[n]) / c->m_j[n];
			q2_z[n] = (d->t2[n] + c->rho_f_j[n]*d->t1[n]) / c->m_j[n];
			v2_z[n] = c->rho_j[n]*d->t1[n] - c->rho_f_j[n]*q2_z[n];
		}
	}
	
    // -------------------------------------------------------------------------
    // \epsilon = \alpha\left(D_x^-v_x + D_z^-v_z\right) + D_x^-q_x + D_z^-q_z
    //
    
    // D_x- v_x
    fftw_execute( d->vx_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->vx_x_i );  // d->t1 now holds D_xv_x
	
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_vx[ipml] = p->b_x[nl]*m[0].dx_vx[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_vx[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_vx[ipml] = p->b_x[nl]*m[0].dx_vx[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_vx[ipml];
        }
    }
    
    // D_z- v_z
    fftw_execute( d->vz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->vz_z_i );  // d->t2 now holds D_zv_z
	
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_vz[ipml] = p->b_z[nl]*m[0].dz_vz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_vz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_vz[ipml] = p->b_z[nl]*m[0].dz_vz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_vz[ipml];
        }
    }
    
    for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] = c->alpha[n] * ( d->t1[n]+d->t2[n] );
    
    // D_x- q_x
    fftw_execute( d->qx_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->qx_x_i );  // d->t3 now holds D_xq_x
	
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_qx[ipml] = p->b_x[nl]*m[0].dx_qx[ipml] + p->c_x[nl]*d->t3[ind];
            d->t3[ind] = d->t3[ind]/p->k_x[nl] + m[2].dx_qx[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_qx[ipml] = p->b_x[nl]*m[0].dx_qx[ipml] + p->c_x[nl]*d->t3[ind];
            d->t3[ind] = d->t3[ind]/p->k_x[nl] + m[2].dx_qx[ipml];
        }
    }

    for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] += d->t3[n];
    
    // D_z- q_z
    fftw_execute( d->qz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->qz_z_i );  // d->t3 now holds D_zq_z
	
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_qz[ipml] = p->b_z[nl]*m[0].dz_qz[ipml] + p->c_z[nl]*d->t3[ind];
            d->t3[ind] = d->t3[ind]/p->k_z[nl] + m[2].dz_qz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_qz[ipml] = p->b_z[nl]*m[0].dz_qz[ipml] + p->c_z[nl]*d->t3[ind];
            d->t3[ind] = d->t3[ind]/p->k_z[nl] + m[2].dz_qz[ipml];
        }
    }
    
    for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] += d->t3[n];
    
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial e}{\partial t} = -\frac{1}{\tau_\sigma}\left[M\left(L+\varphi\right)^{-1}\varphi\epsilon+e\right]
    //
    //      1/tau_s                              is stored into   tau_s
    //      M*varphi/(1+varphi)                  is stored into   varphi
    
    for ( size_t n=0; n<nnodes; ++n ) {
        e2[n] = -c->tau_s[n] * (c->varphi[n]*c->epsilon[n] + e[n]);
        //  store   M*epsilon + e  in  epsilon
        c->epsilon[n] = c->M[n]*c->epsilon[n] + e[n];
    }
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial p}{\partial t} = -\left(M\epsilon+e\right)+s_f
    //
    for ( size_t n=0; n<nnodes; ++n )
        p2[n] = - c->epsilon[n];
	
	
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{xx}}{\partial t} = E D_x^-v_x + (E-2\mu)D_z^-v_z+\alpha(M\epsilon+e)+s_{x}
    //
    // c->E  holds in fact E
    // c->mu holds in fact (E-2mu)
    for ( size_t n=0; n<nnodes; ++n )
        tau2_xx[n] = c->E[n]*d->t1[n] + c->mu[n]*d->t2[n] + c->alpha[n]*c->epsilon[n];
	
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{zz}}{\partial t} = (E-2\mu) D_x^-v_x + E D_z^-v_z+\alpha(M\epsilon+e)+s_{z}
    //
    for ( size_t n=0; n<nnodes; ++n )
        tau2_zz[n] = c->mu[n]*d->t1[n] + c->E[n]*d->t2[n] + c->alpha[n]*c->epsilon[n];
	
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{xz}}{\partial t} = \mu\left(D_z^+v_x + D_x^+v_z\right)+s_{xz}
    //
    //    with \mu interpolated at (i+1/2,j+1/2)
    
    // D_z+ v_x
    fftw_execute( d->vx_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->vx_z_i );  // d->t1 now holds D_zv_x
	
    // D_x+ v_z
    fftw_execute( d->vz_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->vz_x_i );  // d->t2 now holds D_xv_z
	
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_vx[ipml] = p->bh_z[nl]*m[0].dz_vx[ipml] + p->ch_z[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kh_z[nl] + m[2].dz_vx[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_vx[ipml] = p->bH_z[nl]*m[0].dz_vx[ipml] + p->cH_z[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kH_z[nl] + m[2].dz_vx[ipml];
        }
    }    
    
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_vz[ipml] = p->bh_x[nl]*m[0].dx_vz[ipml] + p->ch_x[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kh_x[nl] + m[2].dx_vz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // right
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_vz[ipml] = p->bH_x[nl]*m[0].dx_vz[ipml] + p->cH_x[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kH_x[nl] + m[2].dx_vz[ipml];
        }
    }
    
    for ( size_t n=0; n<nnodes; ++n )
        tau2_xz[n] = c->mu_ij[n] * ( d->t1[n]+d->t2[n] );
    
	// update pml memory variables
	for ( size_t n=0; n<(2*g->ab.np+1)*g->nz2; ++n ) {
		m[2].dx_txx[n] = m[1].dx_txx[n];
		m[2].dx_p[n]   = m[1].dx_p[n];
		m[2].dx_txz[n] = m[1].dx_txz[n];
		m[2].dx_vx[n]  = m[1].dx_vx[n];
		m[2].dx_qx[n]  = m[1].dx_qx[n];
		m[2].dx_vz[n]  = m[1].dx_vz[n];
	}
	for ( size_t n=0; n<(2*g->ab.np+1)*g->nx2; ++n ) {
		m[2].dz_txz[n] = m[1].dz_txz[n];
		m[2].dz_tzz[n] = m[1].dz_tzz[n];
		m[2].dz_p[n]   = m[1].dz_p[n];
		m[2].dz_vz[n]  = m[1].dz_vz[n];
		m[2].dz_qz[n]  = m[1].dz_qz[n];
		m[2].dz_vx[n]  = m[1].dz_vx[n]; 
	}			
}

void propagateVTI(double *W, double *delta, const struct grid *g,
				  const struct computationVariablesVTI *c, struct fftw_data *d) {

    size_t nnodes = g->nx2 * g->nz2;
	
//    double *q_x     = &(W[6*nnodes]);
//    double *q_z     = &(W[7*nnodes]);
    double *e       = &(W[8*nnodes]);
    
    double *tau2_xx = &(delta[0]);
    double *tau2_zz = &(delta[nnodes]);
    double *tau2_xz = &(delta[2*nnodes]);
    double *p2      = &(delta[3*nnodes]);
    double *v2_x    = &(delta[4*nnodes]);
    double *v2_z    = &(delta[5*nnodes]);
    double *q2_x    = &(delta[6*nnodes]);
    double *q2_z    = &(delta[7*nnodes]);
    double *e2      = &(delta[8*nnodes]);
    
	// -------------------------------------------------------------------------
    //
    //  Solve for 
    //
    // D_x^+\tau_{xx} + D_z^-\tau_{xz} = \rho\frac{\partial v_x}{\partial t} + \rho_f\frac{\partial q_x}{\partial t}
    //    and
    // -D_x^+p = \rho_f\frac{\partial v_x}{\partial t} + m_1\frac{\partial q_x}{\partial t} + \frac{\eta}{\kappa_1}q_x
    //
    //  Note:
    //           rho_f^2/rho - m_1                    is stored into   m1
    //           rho_f/rho                            is stored into   rho_f
    //           1/rho                                is stored into   rho
	
    
    // D_x+ tau_xx
    fftw_execute( d->tauxx_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->tauxx_x_i );  // d->t1 now holds D_x tau_xx
    // Dz- tau_xz
    fftw_execute( d->tauxz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->tauxz_z_i );  // d->t2 now holds D_z tau_xz
    
    for ( size_t n=0; n<nnodes; ++n )
        d->t1[n] += d->t2[n];  // d->t1 now holds D_x tau_xx  +  D_z tau_xz
    
    // D_x+ p
    fftw_execute( d->p_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->p_x_i );  // d->t2 now holds D_x p
    
    for ( size_t n=0; n<nnodes; ++n ) {
//        q2_x[n] = (d->t2[n] + c->rho_f_i[n]*d->t1[n] + c->nk1[n]*q_x[n]) / c->m1[n];
        q2_x[n] = (d->t2[n] + c->rho_f_i[n]*d->t1[n]) / c->m1[n];
        v2_x[n] = c->rho_i[n]*d->t1[n] - c->rho_f_i[n]*q2_x[n];
    }    
    
    
    
    // -------------------------------------------------------------------------
    //
    //  Solve for 
    //
    // D_x^+\tau_{xz} + D_z^-\tau_{zz} = \rho\frac{\partial v_z}{\partial t} + \rho_f\frac{\partial q_z}{\partial t}
    //    and
    // -D_z^+p = \rho_f\frac{\partial v_z}{\partial t} + m_3\frac{\partial q_z}{\partial t} + \frac{\eta}{\kappa_3}q_z
	//
    //           rho_f^2/rho - m_3                    is stored into   m3
    
    // D_x- tau_xz
    fftw_execute( d->tauxz_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->tauxz_x_i );  // d->t1 now holds D_x tau_xz
    // D_z+ tau_zz
    fftw_execute( d->tauzz_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->tauzz_z_i );  // d->t2 now holds D_z tau_zz
    
    for ( size_t n=0; n<nnodes; ++n )
        d->t1[n] += d->t2[n];  // d->t1 now holds D_x tau_xz  +  D_z tau_zz
	
    // D_z+ p
    fftw_execute( d->p_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->p_z_i );  // d->t2 now holds D_z p
	
    for ( size_t n=0; n<nnodes; ++n ) {
//        q2_z[n] = (d->t2[n] + c->rho_f_j[n]*d->t1[n] + c->nk3[n]*q_z[n]) / c->m3[n];
        q2_z[n] = (d->t2[n] + c->rho_f_j[n]*d->t1[n]) / c->m3[n];
        v2_z[n] = c->rho_j[n]*d->t1[n] - c->rho_f_j[n]*q2_z[n];
    }    
	
    
	
    // -------------------------------------------------------------------------
    // \epsilon = \alpha_1 D_x^-v_x + \alpha_3 D_z^-v_z + D_x^-q_x + D_z^-q_z
    //
    
    // D_x- v_x
    fftw_execute( d->vx_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->vx_x_i );  // d->t1 now holds D_xv_x
    
    // D_z- v_z
    fftw_execute( d->vz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->vz_z_i );  // d->t2 now holds D_zv_z
    
    for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] = c->alpha1[n] * d->t1[n] + c->alpha3[n] * d->t2[n];
    
    // D_x- q_x
    fftw_execute( d->qx_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->qx_x_i );  // d->t3 now holds D_xq_x
    for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] += d->t3[n];
    
    // D_z- q_z
    fftw_execute( d->qz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->qz_z_i );  // d->t3 now holds D_zq_z
    for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] += d->t3[n];
    
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial e}{\partial t} = -\frac{1}{\tau_\sigma}\left[M\left(L+\varphi\right)^{-1}\varphi\epsilon+e\right]
    //
    //      1/tau_s                              is stored into   tau_s
    //      M*varphi/(1+varphi)                  is stored into   varphi
	
    for ( size_t n=0; n<nnodes; ++n ) {
        e2[n] = -c->tau_s[n] * (c->varphi[n]*c->epsilon[n] + e[n]);
        //  store   M*epsilon + e  in  epsilon
        c->epsilon[n] = c->M[n]*c->epsilon[n] + e[n];
    }
	
    // -------------------------------------------------------------------------
    //
    // \frac{\partial p}{\partial t} = -\left(M\epsilon+e\right)+s_f
    //
    for ( size_t n=0; n<nnodes; ++n )
        p2[n] = - c->epsilon[n];
	
	
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{xx}}{\partial t} = c_{11} D_x^-v_x + c_{13} D_z^-v_z+\alpha_1(M\epsilon+e)+s_{x}
    //
    for ( size_t n=0; n<nnodes; ++n )
        tau2_xx[n] = c->c11[n]*d->t1[n] + c->c13[n]*d->t2[n] + c->alpha1[n]*c->epsilon[n];
	
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{zz}}{\partial t} = c_{13} D_x^-v_x + c_{33} D_z^-v_z+\alpha_3(M\epsilon+e)+s_{z}
    //
    for ( size_t n=0; n<nnodes; ++n )
        tau2_zz[n] = c->c13[n]*d->t1[n] + c->c33[n]*d->t2[n] + c->alpha3[n]*c->epsilon[n];
	
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{xz}}{\partial t} = \mu\left(D_z^+v_x + D_x^+v_z\right)+s_{xz}
    //
    //    with c55 interpolated at (i+1/2,j+1/2)
    
    // D_z+ v_x
    fftw_execute( d->vx_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->vx_z_i );  // d->t1 now holds D_zv_x
	
    // D_x+ v_z
    fftw_execute( d->vz_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->vz_x_i );  // d->t2 now holds D_xv_z
	
    for ( size_t n=0; n<nnodes; ++n )
        tau2_xz[n] = c->c55[n] * ( d->t1[n]+d->t2[n] );
    
}

void propagateVTI_CPML(double *W, double *delta, const struct grid *g,
					   const struct computationVariablesVTI *c,
					   struct fftw_data *d, struct mem_pml *m,
					   struct fac_pml *p) {
	
    size_t nnodes = g->nx2 * g->nz2;
	size_t npml = g->ab.np;

//    double *q_x     = &(W[6*nnodes]);
//    double *q_z     = &(W[7*nnodes]);
    double *e       = &(W[8*nnodes]);
    
    double *tau2_xx = &(delta[0]);
    double *tau2_zz = &(delta[nnodes]);
    double *tau2_xz = &(delta[2*nnodes]);
    double *p2      = &(delta[3*nnodes]);
    double *v2_x    = &(delta[4*nnodes]);
    double *v2_z    = &(delta[5*nnodes]);
    double *q2_x    = &(delta[6*nnodes]);
    double *q2_z    = &(delta[7*nnodes]);
    double *e2      = &(delta[8*nnodes]);
    
	// -------------------------------------------------------------------------
    //
    //  Solve for 
    //
    // D_x^+\tau_{xx} + D_z^-\tau_{xz} = \rho\frac{\partial v_x}{\partial t} + \rho_f\frac{\partial q_x}{\partial t}
    //    and
    // -D_x^+p = \rho_f\frac{\partial v_x}{\partial t} + m_1\frac{\partial q_x}{\partial t} + \frac{\eta}{\kappa_1}q_x
    //
    //  Note:
    //           rho_f^2/rho - m_1                    is stored into   m1
    //           rho_f/rho                            is stored into   rho_f
    //           1/rho                                is stored into   rho
	
    
    // D_x+ tau_xx
    fftw_execute( d->tauxx_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->tauxx_x_i );  // d->t1 now holds D_x tau_xx
    // Dz- tau_xz
    fftw_execute( d->tauxz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->tauxz_z_i );  // d->t2 now holds D_z tau_xz
	
	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j
	size_t ipml = 0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_txx[ipml] = p->bh_x[nl]*m[0].dx_txx[ipml] + p->ch_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kh_x[nl] + m[2].dx_txx[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // right
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_txx[ipml] = p->bH_x[nl]*m[0].dx_txx[ipml] + p->cH_x[nl] *d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kH_x[nl] + m[2].dx_txx[ipml];
        }
    }
    // derivative evaluated at i+1/2,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_txz[ipml] = p->b_z[nl]*m[0].dz_txz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_txz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_txz[ipml] = p->b_z[nl]*m[0].dz_txz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_txz[ipml];
        }
    }
	
    for ( size_t n=0; n<nnodes; ++n )
        d->t1[n] += d->t2[n];  // d->t1 now holds D_x tau_xx  +  D_z tau_xz
    
    // D_x+ p
    fftw_execute( d->p_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->p_x_i );  // d->t2 now holds D_x p
    
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_p[ipml] = p->bh_x[nl]*m[0].dx_p[ipml] + p->ch_x[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kh_x[nl] + m[2].dx_p[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // bottom
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_p[ipml] = p->bH_x[nl]*m[0].dx_p[ipml] + p->cH_x[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kH_x[nl] + m[2].dx_p[ipml];            
        }
    }

	for (size_t i=1; i<g->nx2-1; ++i) {
		for ( size_t j=1; j<g->nz2-1; ++j) {
			size_t n = i*g->nz2+j;
//			q2_x[n] = (d->t2[n] + c->rho_f_i[n]*d->t1[n] + c->nk1[n]*q_x[n]) / c->m1[n];
			q2_x[n] = (d->t2[n] + c->rho_f_i[n]*d->t1[n]) / c->m1[n];
			v2_x[n] = c->rho_i[n]*d->t1[n] - c->rho_f_i[n]*q2_x[n];
		}
	}
    
    
    
    // -------------------------------------------------------------------------
    //
    //  Solve for 
    //
    // D_x^+\tau_{xz} + D_z^-\tau_{zz} = \rho\frac{\partial v_z}{\partial t} + \rho_f\frac{\partial q_z}{\partial t}
    //    and
    // -D_z^+p = \rho_f\frac{\partial v_z}{\partial t} + m_3\frac{\partial q_z}{\partial t} + \frac{\eta}{\kappa_3}q_z
	//
    //           rho_f^2/rho - m_3                    is stored into   m3
    
    // D_x- tau_xz
    fftw_execute( d->tauxz_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->tauxz_x_i );  // d->t1 now holds D_x tau_xz
    // D_z+ tau_zz
    fftw_execute( d->tauzz_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->tauzz_z_i );  // d->t2 now holds D_z tau_zz
    
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_txz[ipml] = p->b_x[nl]*m[0].dx_txz[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_txz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {  // right
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_txz[ipml] = p->b_x[nl]*m[0].dx_txz[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_txz[ipml];
        }
    }
    // derivative evaluated at i,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_tzz[ipml] = p->bh_z[nl]*m[0].dz_tzz[ipml] + p->ch_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kh_z[nl] + m[2].dz_tzz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_tzz[ipml] = p->bH_z[nl]*m[0].dz_tzz[ipml] + p->cH_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kH_z[nl] + m[2].dz_tzz[ipml];
        }
    }
    
    for ( size_t n=0; n<nnodes; ++n )
        d->t1[n] += d->t2[n];  // d->t1 now holds D_x tau_xz  +  D_z tau_zz
	
    // D_z+ p
    fftw_execute( d->p_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->p_z_i );  // d->t2 now holds D_z p
	
    // derivative evaluated at i,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_p[ipml] = p->bh_z[nl]*m[0].dz_p[ipml] + p->ch_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kh_z[nl] + m[2].dz_p[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_p[ipml] = p->bH_z[nl]*m[0].dz_p[ipml] + p->cH_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kH_z[nl] + m[2].dz_p[ipml];
        }
    }
	
	for (size_t i=1; i<g->nx2-1; ++i) {
		for ( size_t j=1; j<g->nz2-1; ++j) {
			size_t n = i*g->nz2+j;
//			q2_z[n] = (d->t2[n] + c->rho_f_j[n]*d->t1[n] + c->nk3[n]*q_z[n]) / c->m3[n];
			q2_z[n] = (d->t2[n] + c->rho_f_j[n]*d->t1[n]) / c->m3[n];
			v2_z[n] = c->rho_j[n]*d->t1[n] - c->rho_f_j[n]*q2_z[n];
		}
	}
	
    
	
    // -------------------------------------------------------------------------
    // \epsilon = \alpha_1 D_x^-v_x + \alpha_3 D_z^-v_z + D_x^-q_x + D_z^-q_z
    //
    
    // D_x- v_x
    fftw_execute( d->vx_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->vx_x_i );  // d->t1 now holds D_xv_x
    
	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_vx[ipml] = p->b_x[nl]*m[0].dx_vx[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_vx[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_vx[ipml] = p->b_x[nl]*m[0].dx_vx[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_vx[ipml];
        }
    }
	
    // D_z- v_z
    fftw_execute( d->vz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->vz_z_i );  // d->t2 now holds D_zv_z

	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_vz[ipml] = p->b_z[nl]*m[0].dz_vz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_vz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_vz[ipml] = p->b_z[nl]*m[0].dz_vz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_vz[ipml];
        }
    }	
    
    for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] = c->alpha1[n] * d->t1[n] + c->alpha3[n] * d->t2[n];
    
    // D_x- q_x
    fftw_execute( d->qx_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->qx_x_i );  // d->t3 now holds D_xq_x
    
	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_qx[ipml] = p->b_x[nl]*m[0].dx_qx[ipml] + p->c_x[nl]*d->t3[ind];
            d->t3[ind] = d->t3[ind]/p->k_x[nl] + m[2].dx_qx[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_qx[ipml] = p->b_x[nl]*m[0].dx_qx[ipml] + p->c_x[nl]*d->t3[ind];
            d->t3[ind] = d->t3[ind]/p->k_x[nl] + m[2].dx_qx[ipml];
        }
    }

	for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] += d->t3[n];
    
    // D_z- q_z
    fftw_execute( d->qz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->qz_z_i );  // d->t3 now holds D_zq_z
	
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_qz[ipml] = p->b_z[nl]*m[0].dz_qz[ipml] + p->c_z[nl]*d->t3[ind];
            d->t3[ind] = d->t3[ind]/p->k_z[nl] + m[2].dz_qz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_qz[ipml] = p->b_z[nl]*m[0].dz_qz[ipml] + p->c_z[nl]*d->t3[ind];
            d->t3[ind] = d->t3[ind]/p->k_z[nl] + m[2].dz_qz[ipml];
        }
    }
    
    for ( size_t n=0; n<nnodes; ++n )
        c->epsilon[n] += d->t3[n];
    
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial e}{\partial t} = -\frac{1}{\tau_\sigma}\left[M\left(L+\varphi\right)^{-1}\varphi\epsilon+e\right]
    //
    //      1/tau_s                              is stored into   tau_s
    //      M*varphi/(1+varphi)                  is stored into   varphi
	
    for ( size_t n=0; n<nnodes; ++n ) {
        e2[n] = -c->tau_s[n] * (c->varphi[n]*c->epsilon[n] + e[n]);
        //  store   M*epsilon + e  in  epsilon
        c->epsilon[n] = c->M[n]*c->epsilon[n] + e[n];
    }
	
    // -------------------------------------------------------------------------
    //
    // \frac{\partial p}{\partial t} = -\left(M\epsilon+e\right)+s_f
    //
    for ( size_t n=0; n<nnodes; ++n )
        p2[n] = - c->epsilon[n];
	
	
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{xx}}{\partial t} = c_{11} D_x^-v_x + c_{13} D_z^-v_z+\alpha_1(M\epsilon+e)+s_{x}
    //
    for ( size_t n=0; n<nnodes; ++n )
        tau2_xx[n] = c->c11[n]*d->t1[n] + c->c13[n]*d->t2[n] + c->alpha1[n]*c->epsilon[n];
	
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{zz}}{\partial t} = c_{13} D_x^-v_x + c_{33} D_z^-v_z+\alpha_3(M\epsilon+e)+s_{z}
    //
    for ( size_t n=0; n<nnodes; ++n )
        tau2_zz[n] = c->c13[n]*d->t1[n] + c->c33[n]*d->t2[n] + c->alpha3[n]*c->epsilon[n];
	
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{xz}}{\partial t} = \mu\left(D_z^+v_x + D_x^+v_z\right)+s_{xz}
    //
    //    with c55 interpolated at (i+1/2,j+1/2)
    
    // D_z+ v_x
    fftw_execute( d->vx_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->vx_z_i );  // d->t1 now holds D_zv_x
	
	// D_x+ v_z
    fftw_execute( d->vz_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->vz_x_i );  // d->t2 now holds D_xv_z
	
	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_vx[ipml] = p->bh_z[nl]*m[0].dz_vx[ipml] + p->ch_z[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kh_z[nl] + m[2].dz_vx[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {

            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_vx[ipml] = p->bH_z[nl]*m[0].dz_vx[ipml] + p->cH_z[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kH_z[nl] + m[2].dz_vx[ipml];
        }
    }    
    
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_vz[ipml] = p->bh_x[nl]*m[0].dx_vz[ipml] + p->ch_x[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kh_x[nl] + m[2].dx_vz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // right
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {

            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_vz[ipml] = p->bH_x[nl]*m[0].dx_vz[ipml] + p->cH_x[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kH_x[nl] + m[2].dx_vz[ipml];
        }
    }
    for ( size_t n=0; n<nnodes; ++n )
        tau2_xz[n] = c->c55[n] * ( d->t1[n]+d->t2[n] );
	
	// update pml memory variables
	for ( size_t n=0; n<(2*g->ab.np+1)*g->nz2; ++n ) {
		m[2].dx_txx[n] = m[1].dx_txx[n];
		m[2].dx_p[n]   = m[1].dx_p[n];
		m[2].dx_txz[n] = m[1].dx_txz[n];
		m[2].dx_vx[n]  = m[1].dx_vx[n];
		m[2].dx_qx[n]  = m[1].dx_qx[n];
		m[2].dx_vz[n]  = m[1].dx_vz[n];
	}
	for ( size_t n=0; n<(2*g->ab.np+1)*g->nx2; ++n ) {
		m[2].dz_txz[n] = m[1].dz_txz[n];
		m[2].dz_tzz[n] = m[1].dz_tzz[n];
		m[2].dz_p[n]   = m[1].dz_p[n];
		m[2].dz_vz[n]  = m[1].dz_vz[n];
		m[2].dz_qz[n]  = m[1].dz_qz[n];
		m[2].dz_vx[n]  = m[1].dz_vx[n]; 
	}			
}

void propagateVE_CPML(double *W, double *delta, const struct grid *g,
					   const struct computationVariablesVE_VTI *c,
					   struct fftw_data *d, struct mem_pml *m,
					   struct fac_pml *p, const int L) {
	
    size_t nnodes = g->nx2 * g->nz2;
	size_t npml = g->ab.np;
    
    double *e1 = &(W[5*nnodes]);
    double *e2 = &(W[(5+L)*nnodes]);
    double *e3 = &(W[(6+L)*nnodes]);
    
    double *tau2_xx = &(delta[0]);
    double *tau2_zz = &(delta[nnodes]);
    double *tau2_xz = &(delta[2*nnodes]);
    double *v2_x    = &(delta[3*nnodes]);
    double *v2_z    = &(delta[4*nnodes]);
    double *e12     = &(delta[5*nnodes]);
    double *e22     = &(delta[(5+L)*nnodes]);
    double *e32     = &(delta[(6+L)*nnodes]);;
    
    
	// -------------------------------------------------------------------------
    //
    //  Solve for
    //
    // D_x^+\tau_{xx} + D_z^-\tau_{xz} = \rho\frac{\partial v_x}{\partial t}
    //
	
    
    // D_x+ tau_xx
    fftw_execute( d->tauxx_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->tauxx_x_i );  // d->t1 now holds D_x tau_xx
    // Dz- tau_xz
    fftw_execute( d->tauxz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->tauxz_z_i );  // d->t2 now holds D_z tau_xz
	
	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j
	size_t ipml = 0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_txx[ipml] = p->bh_x[nl]*m[0].dx_txx[ipml] + p->ch_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kh_x[nl] + m[2].dx_txx[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // right
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_txx[ipml] = p->bH_x[nl]*m[0].dx_txx[ipml] + p->cH_x[nl] *d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kH_x[nl] + m[2].dx_txx[ipml];
        }
    }
    // derivative evaluated at i+1/2,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_txz[ipml] = p->b_z[nl]*m[0].dz_txz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_txz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_txz[ipml] = p->b_z[nl]*m[0].dz_txz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_txz[ipml];
        }
    }
	
    for ( size_t n=0; n<nnodes; ++n )
        d->t1[n] += d->t2[n];  // d->t1 now holds D_x tau_xx  +  D_z tau_xz
    
    
	for (size_t i=1; i<g->nx2-1; ++i) {
		for ( size_t j=1; j<g->nz2-1; ++j) {
			size_t n = i*g->nz2+j;
			v2_x[n] = c->rho_i[n]*d->t1[n];
		}
	}
    
    
    
    // -------------------------------------------------------------------------
    //
    //  Solve for
    //
    // D_x^+\tau_{xz} + D_z^-\tau_{zz} = \rho\frac{\partial v_z}{\partial t}
	//
    
    // D_x- tau_xz
    fftw_execute( d->tauxz_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->tauxz_x_i );  // d->t1 now holds D_x tau_xz
    // D_z+ tau_zz
    fftw_execute( d->tauzz_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->tauzz_z_i );  // d->t2 now holds D_z tau_zz
    
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_txz[ipml] = p->b_x[nl]*m[0].dx_txz[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_txz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {  // right
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_txz[ipml] = p->b_x[nl]*m[0].dx_txz[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_txz[ipml];
        }
    }
    // derivative evaluated at i,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_tzz[ipml] = p->bh_z[nl]*m[0].dz_tzz[ipml] + p->ch_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kh_z[nl] + m[2].dz_tzz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_tzz[ipml] = p->bH_z[nl]*m[0].dz_tzz[ipml] + p->cH_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kH_z[nl] + m[2].dz_tzz[ipml];
        }
    }
    
    for ( size_t n=0; n<nnodes; ++n )
        d->t1[n] += d->t2[n];  // d->t1 now holds D_x tau_xz  +  D_z tau_zz
	
	for (size_t i=1; i<g->nx2-1; ++i) {
		for ( size_t j=1; j<g->nz2-1; ++j) {
			size_t n = i*g->nz2+j;
			v2_z[n] = c->rho_j[n]*d->t1[n];
		}
	}
	
    
    // -------------------------------------------------------------------------
    //
    // \frac{\partial e_{1l}}{\partial t} = \frac{1}{\tau_{\sigma l}^{(1)}}\left[\frac{\left(1+\eta_{1l}^{-1}\right)}{L_1}\left(D_x^-v_x+D_z^-v_z\right)-e_{1l}\right]
    //
	// \frac{\partial e_2}{\partial t} = \frac{1}{2\tau_\sigma^{(2)}}\left[\left(1+\eta_2^{-1}\right)\left(D_x^-v_x+D_z^-v_z\right)-2e_2\right]
	//
    //      1/tau_s                              is stored into   tau_s
    //      (1_eta^-1)/L                         is stored into   eta
	
    // D_x- v_x
    fftw_execute( d->vx_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->vx_x_i );  // d->t1 now holds D_xv_x
    
	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_vx[ipml] = p->b_x[nl]*m[0].dx_vx[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_vx[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_vx[ipml] = p->b_x[nl]*m[0].dx_vx[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_vx[ipml];
        }
    }
	
    // D_z- v_z
    fftw_execute( d->vz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->vz_z_i );  // d->t2 now holds D_zv_z
    
	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_vz[ipml] = p->b_z[nl]*m[0].dz_vz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_vz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_vz[ipml] = p->b_z[nl]*m[0].dz_vz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_vz[ipml];
        }
    }

	for ( size_t m=0; m<L; ++m ) {
		for ( size_t n=0; n<nnodes; ++n ) {
			e12[m*nnodes+n] = c->tau_s_1[m][n] * (c->eta1[m][n]*(d->t1[n]+d->t2[n]) - e1[m*nnodes+n]);
		}
	}
	for ( size_t n=0; n<nnodes; ++n ) {
		e22[n] = 0.5 * c->tau_s_2[n] * (c->eta2[n]*(d->t1[n]+d->t2[n]) - 2.*e2[n]);
	}
	
	
	
    // -------------------------------------------------------------------------
    //
    // \frac{\partial \tau_{xx}}{\partial t} = c_{11} D_x^-v_x + c_{13} D_z^-v_z + K^0\epsilon_1 + 2c_{55}^0\epsilon_2
    //
    // \frac{\partial \tau_{zz}}{\partial t} = c_{13} D_x^-v_x + c_{33} D_z^-v_z + K^0\epsilon_1 - 2c_{55}^0\epsilon_2
    //
    for ( size_t n=0; n<nnodes; ++n ) {
		double sum_e1 = 0.0;
		for ( size_t m=0; m<L; ++m ) {
			sum_e1 += e1[m+nnodes+n];
		}
        tau2_xx[n] = c->c11[n]*d->t1[n] + c->c13[n]*d->t2[n] + c->K0[n]*sum_e1 + 2.*c->c55_0[n]*e2[n];
        tau2_zz[n] = c->c13[n]*d->t1[n] + c->c33[n]*d->t2[n] + c->K0[n]*sum_e1 - 2.*c->c55_0[n]*e2[n];
	}
    
    
    // -------------------------------------------------------------------------
    //
	// \frac{\partial e_3}{\partial t} = \frac{1}{\tau_\sigma^{(2)}}\left[\left(1+\eta_2^{-1}\right)\left(D_z^+v_x+D_x^+v_z\right)-e_3\right]
	//
    // \frac{\partial \tau_{xz}}{\partial t} = c_{55}\left(D_z^+v_x + D_x^+v_z\right) + c_{55}^0\epsilon_3
    //
	//    with c55 interpolated at (i+1/2,j+1/2)

    // D_z+ v_x
    fftw_execute( d->vx_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->vx_z_i );  // d->t1 now holds D_zv_x
	
	// D_x+ v_z
    fftw_execute( d->vz_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->vz_x_i );  // d->t2 now holds D_xv_z
	
	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_vx[ipml] = p->bh_z[nl]*m[0].dz_vx[ipml] + p->ch_z[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kh_z[nl] + m[2].dz_vx[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_vx[ipml] = p->bH_z[nl]*m[0].dz_vx[ipml] + p->cH_z[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kH_z[nl] + m[2].dz_vx[ipml];
        }
    }
    
    // update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_vz[ipml] = p->bh_x[nl]*m[0].dx_vz[ipml] + p->ch_x[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kh_x[nl] + m[2].dx_vz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // right
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_vz[ipml] = p->bH_x[nl]*m[0].dx_vz[ipml] + p->cH_x[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->kH_x[nl] + m[2].dx_vz[ipml];
        }
    }
	
	for ( size_t n=0; n<nnodes; ++n ) {
		e32[n] = c->tau_s_2_ij[n] * (c->eta2_ij[n]*(d->t1[n]+d->t2[n]) - e3[n]);
	}
	
    for ( size_t n=0; n<nnodes; ++n )
        tau2_xz[n] = c->c55[n] * ( d->t1[n]+d->t2[n] ) + c->c55_0_ij[n]*e3[n];
	
	// update pml memory variables
	for ( size_t n=0; n<(2*g->ab.np+1)*g->nz2; ++n ) {
		m[2].dx_txx[n] = m[1].dx_txx[n];
		m[2].dx_txz[n] = m[1].dx_txz[n];
		m[2].dx_vx[n]  = m[1].dx_vx[n];
		m[2].dx_vz[n]  = m[1].dx_vz[n];
	}
	for ( size_t n=0; n<(2*g->ab.np+1)*g->nx2; ++n ) {
		m[2].dz_txz[n] = m[1].dz_txz[n];
		m[2].dz_tzz[n] = m[1].dz_tzz[n];
		m[2].dz_vz[n]  = m[1].dz_vz[n];
		m[2].dz_vx[n]  = m[1].dz_vx[n];
	}			
}

void propagateVE_SH_CPML(double *W, double *delta, const struct grid *g,
						 const struct computationVariablesVE_SH_VTI *c,
						 struct fftw_data_ve_sh *d, struct mem_pml_sh *m,
						 struct fac_pml *p) {
	
    size_t nnodes = g->nx2 * g->nz2;
	size_t npml = g->ab.np;
    
    double *e1 = &(W[3*nnodes]);
    double *e2 = &(W[4*nnodes]);
    
    double *tau2_xy = &(delta[0]);
    double *tau2_yz = &(delta[nnodes]);
    double *v2_y    = &(delta[2*nnodes]);
    double *e12     = &(delta[3*nnodes]);
    double *e22     = &(delta[4*nnodes]);
    
    
	// -------------------------------------------------------------------------
    //
    //  Solve for
    //
    // D_x^-\tau_{xy} + D_z^-\tau_{yz} = \rho\frac{\partial v_y}{\partial t}
    //
	
    
    // D_x+ tau_xx
    fftw_execute( d->tauxy_x_f );
    partial_x_backward_sh(d, g);
    fftw_execute( d->tauxy_x_i );  // d->t1 now holds D_x tau_xy
    // Dz- tau_xz
    fftw_execute( d->tauyz_z_f );
    partial_z_backward_sh(d, g);
    fftw_execute( d->tauyz_z_i );  // d->t2 now holds D_z tau_yz
	
	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j
	size_t ipml = 0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // left
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = nl*g->nz2 + j;
            m[1].dx_txy[ipml] = p->b_x[nl]*m[0].dx_txy[ipml] + p->c_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_txy[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // right
        for ( size_t j=0; j<g->nz2; ++j, ++ipml ) {
            
            size_t ind = (g->nx2-1-nl)*g->nz2 + j;
            m[1].dx_txy[ipml] = p->b_x[nl]*m[0].dx_txy[ipml] + p->c_x[nl] *d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->k_x[nl] + m[2].dx_txy[ipml];
        }
    }
    // derivative evaluated at i,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_tyz[ipml] = p->b_z[nl]*m[0].dz_tyz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_tyz[ipml];
		}
	}
	for ( size_t nl=0; nl<npml; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_tyz[ipml] = p->b_z[nl]*m[0].dz_tyz[ipml] + p->c_z[nl]*d->t2[ind];
            d->t2[ind] = d->t2[ind]/p->k_z[nl] + m[2].dz_tyz[ipml];
        }
    }
	
    for ( size_t n=0; n<nnodes; ++n )
        d->t1[n] += d->t2[n];  // d->t1 now holds D_x tau_xy  +  D_z tau_yz
    
    
	for (size_t i=1; i<g->nx2-1; ++i) {
		for ( size_t j=1; j<g->nz2-1; ++j) {
			size_t n = i*g->nz2+j;
			v2_y[n] = c->rho[n]*d->t1[n];
		}
	}
    
    
    
    // -------------------------------------------------------------------------
    //
	// \frac{\partial e_1}{\partial t} = \frac{1}{\tau_\sigma^{(1)}}\left[\left(1+\eta_1^{-1}\right)\left(D_z^+v_y\right)-e_1\right]
	//
    // \frac{\partial \tau_{yz}}{\partial t} = c_{44}\left(D_z^+v_y\right) + c_{44}^0\epsilon_1
    //
	//    with c44 interpolated at (i,j+1/2)
	
    // D_z+ v_y
    fftw_execute( d->vy_z_f );
    partial_z_forward_sh(d, g);
    fftw_execute( d->vy_z_i );  // d->t1 now holds D_zv_y
	
	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i,j+1/2
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dz_vy[ipml] = p->bh_z[nl]*m[0].dz_vy[ipml] + p->ch_z[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kh_z[nl] + m[2].dz_vy[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dz_vy[ipml] = p->bH_z[nl]*m[0].dz_vy[ipml] + p->cH_z[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kH_z[nl] + m[2].dz_vy[ipml];
        }
    }
    
	
	for ( size_t n=0; n<nnodes; ++n ) {
		e12[n] = c->tau_s_1[n] * (c->eta1[n]*(d->t1[n]) - e1[n]);
	}
	
    for ( size_t n=0; n<nnodes; ++n )
        tau2_yz[n] = c->c44[n] * ( d->t1[n] ) + c->c44_0[n]*e1[n];
	
		
	
    // -------------------------------------------------------------------------
    //
	// \frac{\partial e_2}{\partial t} = \frac{1}{\tau_\sigma^{(2)}}\left[\left(1+\eta_2^{-1}\right)\left(D_x^+v_y\right)-e_2\right]
	//
    // \frac{\partial \tau_{xy}}{\partial t} = c_{66}\left(D_x^+v_y\right) + c_{66}^0\epsilon_2
    //
	//    with c66 interpolated at (i+1/2,j)
	
    // D_x+ v_x
    fftw_execute( d->vy_x_f );
    partial_x_forward_sh(d, g);
    fftw_execute( d->vy_x_i );  // d->t1 now holds D_xv_y
	
	// update memory variables and replace spatial derivative in PML region
    // derivative evaluated at i+1/2,j
	ipml=0;
    for ( size_t nl=0; nl<npml; ++nl ) {  // top
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2+nl;
            m[1].dx_vy[ipml] = p->bh_x[nl]*m[0].dx_vy[ipml] + p->ch_x[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kh_x[nl] + m[2].dx_vy[ipml];
		}
	}
	for ( size_t nl=0; nl<npml+1; ++nl ) {  // bottom
        for ( size_t i=0; i<g->nx2; ++i, ++ipml ) {
            
            size_t ind = i*g->nz2 + g->nz2-1-nl;
            m[1].dx_vy[ipml] = p->bH_z[nl]*m[0].dx_vy[ipml] + p->cH_z[nl]*d->t1[ind];
            d->t1[ind] = d->t1[ind]/p->kH_z[nl] + m[2].dx_vy[ipml];
        }
    }
    
	
	for ( size_t n=0; n<nnodes; ++n ) {
		e22[n] = c->tau_s_2[n] * (c->eta2[n]*(d->t1[n]) - e2[n]);
	}
	
    for ( size_t n=0; n<nnodes; ++n )
        tau2_xy[n] = c->c66[n] * ( d->t1[n] ) + c->c66_0[n]*e2[n];
	
	
		
	// update pml memory variables
	for ( size_t n=0; n<(2*g->ab.np+1)*g->nz2; ++n ) {
		m[2].dx_txy[n] = m[1].dx_txy[n];
		m[2].dx_vy[n]  = m[1].dx_vy[n];
	}
	for ( size_t n=0; n<(2*g->ab.np+1)*g->nx2; ++n ) {
		m[2].dz_tyz[n] = m[1].dz_tyz[n];
		m[2].dz_vy[n]  = m[1].dz_vy[n];
	}
}

void partial_x_forward(struct fftw_data *d, const struct grid *g) {
    for ( size_t n1=0, l=0; n1<g->nz2; ++n1 )
        for ( size_t n2=0; n2<g->nx2/2+1; ++n2, ++l )
            d->o_x[l] *= d->kx_f[n2];    
}

void partial_x_backward(struct fftw_data *d, const struct grid *g) {
    for ( size_t n1=0, l=0; n1<g->nz2; ++n1 )
        for ( size_t n2=0; n2<g->nx2/2+1; ++n2, ++l )
            d->o_x[l] *= d->kx_b[n2];    
}

void partial_z_forward(struct fftw_data *d, const struct grid *g) {
    for ( size_t n1=0, l=0; n1<g->nx2; ++n1 )
        for ( size_t n2=0; n2<g->nz2/2+1; ++n2, ++l )
            d->o_z[l] *= d->kz_f[n2];
}    

void partial_z_backward(struct fftw_data *d, const struct grid *g) {
    for ( size_t n1=0, l=0; n1<g->nx2; ++n1 )
        for ( size_t n2=0; n2<g->nz2/2+1; ++n2, ++l )
            d->o_z[l] *= d->kz_b[n2];
}

void partial_x_forward_sh(struct fftw_data_ve_sh *d, const struct grid *g) {
    for ( size_t n1=0, l=0; n1<g->nz2; ++n1 )
        for ( size_t n2=0; n2<g->nx2/2+1; ++n2, ++l )
            d->o_x[l] *= d->kx_f[n2];
}

void partial_x_backward_sh(struct fftw_data_ve_sh *d, const struct grid *g) {
    for ( size_t n1=0, l=0; n1<g->nz2; ++n1 )
        for ( size_t n2=0; n2<g->nx2/2+1; ++n2, ++l )
            d->o_x[l] *= d->kx_b[n2];
}

void partial_z_forward_sh(struct fftw_data_ve_sh *d, const struct grid *g) {
    for ( size_t n1=0, l=0; n1<g->nx2; ++n1 )
        for ( size_t n2=0; n2<g->nz2/2+1; ++n2, ++l )
            d->o_z[l] *= d->kz_f[n2];
}

void partial_z_backward_sh(struct fftw_data_ve_sh *d, const struct grid *g) {
    for ( size_t n1=0, l=0; n1<g->nx2; ++n1 )
        for ( size_t n2=0; n2<g->nz2/2+1; ++n2, ++l )
            d->o_z[l] *= d->kz_b[n2];
}

void update_vx(double *vx, const double *txx, const double *txz, const double *b,
               const double *txx_x, const double *txx_z,
			   const double *txz_x, const double *txz_z,
               const struct grid *g, const double dt, const double h,
               const size_t Npml) {
    
    register double txx_pml, txz_pml, dt1, dt2;
    size_t off2 = Npml*Npml;
    txx_pml = txx_x[off2 + (Npml-1)*g->nz] + txx_z[off2 + (Npml-1)*g->nz];
    off2 = 2*Npml*Npml + Npml*g->nz;
    txz_pml = txz_x[off2 + Npml-1] + txz_z[off2 + Npml-1];
    dt1 = (txx[0]-txx_pml) / h;
    dt2 = (txz[0]-txz_pml) / h;
    
    vx[0] += dt*b[0] * (dt1 + dt2);
    
    off2 = Npml*Npml;
    for (size_t j=1; j<g->nz; ++j) {
        txx_pml = txx_x[off2 + (Npml-1)*g->nz+j] + txx_z[off2 + (Npml-1)*g->nz+j];
        dt1 = (txx[j]-txx_pml) / h;
        dt2 = (txz[j]-txz[j-1]  ) / h;
        
        vx[j] += dt*b[j] * (dt1 + dt2);
    }
    off2 = 2*Npml*Npml + Npml*g->nz;
    for (size_t i=1; i<g->nx; ++i) {
        txz_pml = txz_x[off2 + i*Npml+Npml-1] + txz_z[off2 + i*Npml+Npml-1];
        dt1 = (txx[i*g->nz]-txx[(i-1)*g->nz]) / h;
        dt2 = (txz[i*g->nz]-txz_pml  ) / h;
        
        vx[i*g->nz] += dt*b[i*g->nz] * (dt1 + dt2);
        
        for (size_t j=1; j<g->nz; ++j) {
            dt1 = (txx[i*g->nz+j]-txx[(i-1)*g->nz+j]) / h;
            dt2 = (txz[i*g->nz+j]-txz[i*g->nz+j-1]  ) / h;
            
            vx[i*g->nz+j] += dt*b[i*g->nz+j] * (dt1 + dt2);
        }
    }
}

void update_vz(double *vz, const double *txz, const double *tzz, const double *bm,
               const double *txz_x, const double *txz_z,
			   const double *tzz_x, const double *tzz_z,
               const struct grid *g, const double dt, const double h,
               const size_t Npml) {
	
    register double dt1, dt2, txz_pml, tzz_pml;
    size_t off2 = 2*Npml*Npml + Npml*g->nz + Npml*g->nx;
    for (size_t i=0; i<g->nx-1; ++i) {
        for (size_t j=0; j<g->nz-1; ++j) {
            dt1 = (txz[(i+1)*g->nz+j] - txz[i*g->nz+j]) / h;
            dt2 = (tzz[i*g->nz+j+1]   - tzz[i*g->nz+j]) / h;
            
            vz[i*g->nz+j] += dt*bm[i*g->nz+j] * (dt1 + dt2);
        }
        tzz_pml = tzz_x[off2 + i*Npml] + tzz_z[off2 + i*Npml];
        dt1 = (txz[(i+1)*g->nz+g->nz-1] - txz[i*g->nz+g->nz-1]) / h;
        dt2 = (tzz_pml                  - tzz[i*g->nz+g->nz-1]) / h;
        
        vz[i*g->nz+g->nz-1] += dt*bm[i*g->nz+g->nz-1] * (dt1 + dt2);
    }
    off2 = 3*Npml*Npml + Npml*g->nz + 2*Npml*g->nx;
    for (size_t j=0; j<g->nz-1; ++j) { // i = g->nx-1
        txz_pml = txz_x[off2 + j] + txz_z[off2 + j];
        dt1 = (txz_pml                - txz[(g->nx-1)*g->nz+j]) / h;
        dt2 = (tzz[(g->nx-1)*g->nz+j+1] - tzz[(g->nx-1)*g->nz+j]) / h;
        
        vz[(g->nx-1)*g->nz+j] += dt*bm[(g->nx-1)*g->nz+j] * (dt1 + dt2);
    }
    txz_pml = txz_x[off2 + g->nz-1] + txz_z[off2 + g->nz-1];
    off2 = 2*Npml*Npml + Npml*g->nz + Npml*g->nx;
    tzz_pml = tzz_x[off2 + (g->nx-1)*Npml] + tzz_z[off2 + (g->nx-1)*Npml];
    dt1 = (txz_pml - txz[g->nx*g->nz-1]) / h;
    dt2 = (tzz_pml - tzz[g->nx*g->nz-1]) / h;
    
    vz[g->nx*g->nz-1] += dt*bm[g->nx*g->nz-1] * (dt1 + dt2);
}

void update_txxzz(double *txx, double *tzz, const double *vx, const double *vz,
                  const double *la, const double *l2m,
                  const double *vx_x, const double *vx_z, const double *vz_x, const double *vz_z,
                  const struct grid *g, const double dt, const double h,
                  const size_t Npml) {
    
    register double vx_pml, vz_pml, dvx_x, dvz_z;
    size_t off2 = 2*Npml*Npml + Npml*g->nz;
    for (size_t i=0; i<g->nx-1; ++i) {
        vz_pml = vz_x[off2 + i*Npml+Npml-1] + vz_z[off2 + i*Npml+Npml-1];
        dvx_x = dt*(vx[(i+1)*g->nz]  - vx[i*g->nz]) / h;
        dvz_z = dt*(vz[i*g->nz]      - vz_pml) / h;
        
        txx[i*g->nz] += l2m[i*g->nz]*dvx_x + la[i*g->nz]*dvz_z;
        tzz[i*g->nz] += l2m[i*g->nz]*dvz_z + la[i*g->nz]*dvx_x;
        for (size_t j=1; j<g->nz; ++j) {
            dvx_x = dt*(vx[(i+1)*g->nz+j]-vx[i*g->nz+j])   / h;
            dvz_z = dt*(vz[i*g->nz+j]    -vz[i*g->nz+j-1]) / h;
            
            txx[i*g->nz+j] += l2m[i*g->nz+j]*dvx_x + la[i*g->nz+j]*dvz_z;
            tzz[i*g->nz+j] += l2m[i*g->nz+j]*dvz_z + la[i*g->nz+j]*dvx_x;
        }
    }
    off2 = 3*Npml*Npml + Npml*g->nz + 2*Npml*g->nx;
    vx_pml = vx_x[off2] + vx_z[off2];
    dvx_x = dt*(vx_pml        - vx[(g->nx-1)*g->nz])   / h;
    off2 = 2*Npml*Npml + Npml*g->nz;
    vz_pml = vz_x[off2 + (g->nx-1)*Npml+Npml-1] + vz_z[off2 + (g->nx-1)*Npml+Npml-1];
    dvz_z = dt*(vz[(g->nx-1)*g->nz] - vz_pml) / h;
    
    txx[(g->nx-1)*g->nz] += l2m[(g->nx-1)*g->nz]*dvx_x + la[(g->nx-1)*g->nz]*dvz_z;
    tzz[(g->nx-1)*g->nz] += l2m[(g->nx-1)*g->nz]*dvz_z + la[(g->nx-1)*g->nz]*dvx_x;
    
    off2 = 3*Npml*Npml + Npml*g->nz + 2*Npml*g->nx;
    for (size_t j=1; j<g->nz; ++j) {
        vx_pml = vx_x[off2+j] + vx_z[off2+j];
        dvx_x = dt*(vx_pml              - vx[(g->nx-1)*g->nz+j])   / h;
        dvz_z = dt*(vz[(g->nx-1)*g->nz+j] - vz[(g->nx-1)*g->nz+j-1]) / h;
        
        txx[(g->nx-1)*g->nz+j] += l2m[(g->nx-1)*g->nz+j]*dvx_x + la[(g->nx-1)*g->nz+j]*dvz_z;
        tzz[(g->nx-1)*g->nz+j] += l2m[(g->nx-1)*g->nz+j]*dvz_z + la[(g->nx-1)*g->nz+j]*dvx_x;
    }
}

void update_txz(double *txz, const double *vx, const double *vz, const double *mu,
                const double *vx_x, const double *vx_z, const double *vz_x, const double *vz_z, 
                const struct grid *g, const double dt, const double h,
                const size_t Npml) {

    register double vx_pml, vz_pml, dvx_z, dvz_x;
    size_t off2 = 2*Npml*Npml + Npml*g->nz + Npml*g->nx;
    for (size_t i=1; i<g->nx; ++i) {
        for (size_t j=0; j<g->nz-1; ++j) {
            dvx_z = dt*(vx[i*g->nz+j+1] - vx[i*g->nz+j]) / h;
            dvz_x = dt*(vz[i*g->nz+j]   - vz[(i-1)*g->nz+j]) / h;
            
            txz[i*g->nz+j] += mu[i*g->nz+j] * (dvx_z + dvz_x);
        }
        vx_pml = vx_x[off2 + i*Npml] + vx_z[off2 + i*Npml];
        dvx_z = dt*(vx_pml            - vx[i*g->nz+g->nz-1]) / h;
        dvz_x = dt*(vz[i*g->nz+g->nz-1] - vz[(i-1)*g->nz+g->nz-1]) / h;
        
        txz[i*g->nz+g->nz-1] += mu[i*g->nz+g->nz-1] * (dvx_z + dvz_x);
    }
    off2 = 2*Npml*Npml + Npml*g->nz + Npml*g->nx;
    vx_pml = vx_x[off2] + vx_z[off2];
    dvx_z = dt*(vx_pml   - vx[g->nz-1]) / h;
    off2 = Npml*Npml;
    vz_pml = vz_x[off2 + (Npml-1)*g->nz+g->nz-1] + vz_z[off2 + (Npml-1)*g->nz+g->nz-1];
    dvz_x = dt*(vz[g->nz-1] - vz_pml) / h;
    
    txz[g->nz-1] += mu[g->nz-1] * (dvx_z + dvz_x);
    
    for (size_t j=0; j<g->nz-1; ++j) {
        vz_pml = vz_x[off2 + (Npml-1)*g->nz+j] + vz_z[off2 + (Npml-1)*g->nz+j];
        dvx_z = dt*(vx[j+1] - vx[j]) / h;
        dvz_x = dt*(vz[j]   - vz_pml) / h;
        
        txz[j] += mu[j] * (dvx_z + dvz_x);
    }
}

void compute_div(double *div, const struct grid *g, struct fftw_data *d) {
    // D_x- v_x
    fftw_execute( d->vx_x_f );
    partial_x_backward(d, g);
    fftw_execute( d->vx_x_i );  // d->t1 now holds D_xv_x

    // D_z- v_z
    fftw_execute( d->vz_z_f );
    partial_z_backward(d, g);
    fftw_execute( d->vz_z_i );  // d->t2 now holds D_zv_z

    size_t nnodes = g->nx2 * g->nz2;
    for ( size_t n=0; n<nnodes; ++n )
        div[n] = d->t1[n] + d->t2[n];

}
void compute_curl(double *curl, const struct grid *g, struct fftw_data *d) {
    // D_z+ v_x
    fftw_execute( d->vx_z_f );
    partial_z_forward(d, g);
    fftw_execute( d->vx_z_i );  // d->t1 now holds D_zv_x
	
    // D_x+ v_z
    fftw_execute( d->vz_x_f );
    partial_x_forward(d, g);
    fftw_execute( d->vz_x_i );  // d->t2 now holds D_xv_z

    size_t nnodes = g->nx2 * g->nz2;
    for ( size_t n=0; n<nnodes; ++n )
        curl[n] = d->t2[n] - d->t1[n];
}

