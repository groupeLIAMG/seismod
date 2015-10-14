/*
 *  bc.c
 *
 *  Created by Bernard Giroux on 10-11-27.
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
#include <stdlib.h>
#include "bc.h"

void fill_strips(const struct grid *g, double *dvx, double *dvz,
                 double *dqx, double *dqz, const double *W) {
    
    size_t nnodes = (g->nx + 2*g->bc.np)*(g->nz + 2*g->bc.np);
    
    const double *v_x    = &(W[4*nnodes]);
    const double *v_z    = &(W[5*nnodes]);
    const double *q_x    = &(W[6*nnodes]);
    const double *q_z    = &(W[7*nnodes]);
    
    size_t id = 0;
    // uppermost band, from top to bottom
    for ( size_t j=0; j<g->bc.np; ++j ) {
        for ( size_t i=j; i<g->nx2-j; ++i, ++id ) {

            size_t ij = i*g->nz2+j;
            dvx[ id ] = v_x[ ij ] ;
            dvz[ id ] = v_z[ ij ] ;
            dqx[ id ] = q_x[ ij ] ;
            dqz[ id ] = q_z[ ij ] ;
        }
    }
    // lowermost band, from bottom to top
    for ( size_t j=0; j<g->bc.np; ++j ) {
        for ( size_t i=j; i<g->nx2-j; ++i, ++id ) {
            
            size_t ij = i*g->nz2+(g->nz2-j-1);
            dvx[ id ] = v_x[ ij ] ;
            dvz[ id ] = v_z[ ij ] ;
            dqx[ id ] = q_x[ ij ] ;
            dqz[ id ] = q_z[ ij ] ;
        }
    }
    // leftmost band, from left to right
    for ( size_t i=0; i<g->bc.np; ++i ) {
        for ( size_t j=i+1; j<g->nz2-j-2; ++j, ++id ) {
            
            size_t ij = i*g->nz2+j;
            dvx[ id ] = v_x[ ij ] ;
            dvz[ id ] = v_z[ ij ] ;
            dqx[ id ] = q_x[ ij ] ;
            dqz[ id ] = q_z[ ij ] ;
        }
    }
    // rightmost band, from right to left
    for ( size_t i=0; i<g->bc.np; ++i ) {
        for ( size_t j=i+1; j<g->nz2-j-2; ++j, ++id ) {

            size_t ij = (g->nx2-i-1)*g->nz2+j;
            dvx[ id ] = v_x[ ij ] ;
            dvz[ id ] = v_z[ ij ] ;
            dqx[ id ] = q_x[ ij ] ;
            dqz[ id ] = q_z[ ij ] ;
        }
    }
}


void gb( struct grid *g ) {

    g->bc.gobx = (double *) malloc( g->bc.np*sizeof(double) );
    g->bc.gobz = (double *) malloc( g->bc.np*sizeof(double) );
    
    double U0 = 1.4*g->bc.vmax/g->dx;
    for (size_t i=0; i<g->bc.np; ++i) {
        double s = i*g->bc.alfa;
        double fac = 0.5 * ( exp(s)+exp(-s) );
        g->bc.gobx[i] = U0 / (fac*fac);
    }
    U0 = 1.4*g->bc.vmax/g->dz;
    for (size_t i=0; i<g->bc.np; ++i) {
        double s = i*g->bc.alfa;
        double fac = 0.5 * ( exp(s)+exp(-s) );
        g->bc.gobz[i] = U0 / (fac*fac);
    }
}

void gb2( struct grid *g, const double F ) {
    
    g->bc.gobx = (double *) malloc( g->bc.np*sizeof(double) );
    g->bc.gobz = (double *) malloc( g->bc.np*sizeof(double) );
    
	// Cerjan et al., 1985
	double a = sqrt(-log(F))/(g->bc.np);
    for (size_t i=0; i<g->bc.np; ++i) {
        g->bc.gobz[i] = g->bc.gobx[i] = exp(-(a*(g->bc.np-i))*(a*(g->bc.np-i)));
    }
}

void gb3( struct grid *g ) {
    
    g->bc.gobx = (double *) malloc( g->bc.np*sizeof(double) );
    g->bc.gobz = (double *) malloc( g->bc.np*sizeof(double) );
    
    double d = 4.5 / (g->bc.np-1.);
    for (size_t i=0; i<g->bc.np; ++i) {
        g->bc.gobz[i] = g->bc.gobx[i] = 0.5 + 0.5*erf(d*i - 2.);
    }
}

void tbc(const struct grid *g, const double *dvx, const double *dvz,
         const double *dqx, const double *dqz, double *delta) {

    size_t nnodes = (g->nx + 2*g->bc.np)*(g->nz + 2*g->bc.np);
    
    double *v_x    = &(delta[4*nnodes]);
    double *v_z    = &(delta[5*nnodes]);
    double *q_x    = &(delta[6*nnodes]);
    double *q_z    = &(delta[7*nnodes]);
    
    size_t id = 0;
    // uppermost band, from top to bottom
    for ( size_t j=0; j<g->bc.np; ++j ) {
        for ( size_t i=j; i<g->nx2-j; ++i, ++id ) {
            
            size_t ij = i*g->nz2+j;
            v_x[ ij ] -= g->bc.gobz[j]*dvx[ id ];
            v_z[ ij ] -= g->bc.gobz[j]*dvz[ id ];
            q_x[ ij ] -= g->bc.gobz[j]*dqx[ id ];
            q_z[ ij ] -= g->bc.gobz[j]*dqz[ id ];
        }
    }
    // lowermost band, from bottom to top
    for ( size_t j=0; j<g->bc.np; ++j ) {
        for ( size_t i=j; i<g->nx2-j; ++i, ++id ) {
            
            size_t ij = i*g->nz2+(g->nz2-j-1);
            v_x[ ij ] -= g->bc.gobz[j]*dvx[ id ];
            v_z[ ij ] -= g->bc.gobz[j]*dvz[ id ];
            q_x[ ij ] -= g->bc.gobz[j]*dqx[ id ];
            q_z[ ij ] -= g->bc.gobz[j]*dqz[ id ];
        }
    }
    // leftmost band, from left to right
    for ( size_t i=0; i<g->bc.np; ++i ) {
        for ( size_t j=i+1; j<g->nz2-i-1; ++j, ++id ) {
            
            size_t ij = i*g->nz2+j;
            v_x[ ij ] -= g->bc.gobx[i]*dvx[ id ];
            v_z[ ij ] -= g->bc.gobx[i]*dvz[ id ];
            q_x[ ij ] -= g->bc.gobx[i]*dqx[ id ];
            q_z[ ij ] -= g->bc.gobx[i]*dqz[ id ];
        }
    }
    // rightmost band, from right to left
    for ( size_t i=0; i<g->bc.np; ++i ) {
        for ( size_t j=i+1; j<g->nz2-i-1; ++j, ++id ) {
            
            size_t ij = (g->nx2-i-1)*g->nz2+j;
            v_x[ ij ] -= g->bc.gobx[i]*dvx[ id ];
            v_z[ ij ] -= g->bc.gobx[i]*dvz[ id ];
            q_x[ ij ] -= g->bc.gobx[i]*dqx[ id ];
            q_z[ ij ] -= g->bc.gobx[i]*dqz[ id ];
        }
    }    
}


void tbc2(const struct grid *g, double *delta) {
    
    size_t nnodes = (g->nx + 2*g->bc.np)*(g->nz + 2*g->bc.np);
    
    double *tau_xx = &(delta[0]);
    double *tau_zz = &(delta[nnodes]);
    double *tau_xz = &(delta[2*nnodes]);
    double *p      = &(delta[3*nnodes]);
    double *v_x    = &(delta[4*nnodes]);
    double *v_z    = &(delta[5*nnodes]);
    double *q_x    = &(delta[6*nnodes]);
    double *q_z    = &(delta[7*nnodes]);
    
//	FILE *fid = fopen("tbc.dat","wt");
	
    // uppermost band, from top to bottom
    for ( size_t j=0; j<g->bc.np; ++j ) {
        for ( size_t i=j; i<g->nx2-j; ++i ) {
//            fprintf(fid,"%zd  %zd  %lf\n", i, j, g->bc.gobz[j]);
			
            size_t ij = i*g->nz2+j;
//            tau_xx[ ij ] *= g->bc.gobz[j];
//            tau_xx[ ij ] *= g->bc.gobz[j];
//            tau_xz[ ij ] *= g->bc.gobz[j];
//            p[ ij ]      *= g->bc.gobz[j];
            v_x[ ij ]    *= g->bc.gobz[j];
            v_z[ ij ]    *= g->bc.gobz[j];
            q_x[ ij ]    *= g->bc.gobz[j];
            q_z[ ij ]    *= g->bc.gobz[j];
        }
    }
    // lowermost band, from bottom to top
    for ( size_t j=0; j<g->bc.np; ++j ) {
        for ( size_t i=j; i<g->nx2-j; ++i ) {
//            fprintf(fid,"%zd  %zd  %lf\n", i, g->nz2-j-1, g->bc.gobz[j]);
            
            size_t ij = i*g->nz2+(g->nz2-j-1);
//            tau_xx[ ij ] *= g->bc.gobz[j];
//            tau_xx[ ij ] *= g->bc.gobz[j];
//            tau_xz[ ij ] *= g->bc.gobz[j];
//            p[ ij ]      *= g->bc.gobz[j];
            v_x[ ij ]    *= g->bc.gobz[j];
            v_z[ ij ]    *= g->bc.gobz[j];
            q_x[ ij ]    *= g->bc.gobz[j];
            q_z[ ij ]    *= g->bc.gobz[j];
        }
    }
    // leftmost band, from left to right
    for ( size_t i=0; i<g->bc.np; ++i ) {
        for ( size_t j=i+1; j<g->nz2-i-1; ++j ) {
//            fprintf(fid,"%zd  %zd  %lf\n", i, j, g->bc.gobx[i]);
			
            size_t ij = i*g->nz2+j;
//            tau_xx[ ij ] *= g->bc.gobx[i];
//            tau_xx[ ij ] *= g->bc.gobx[i];
//            tau_xz[ ij ] *= g->bc.gobx[i];
//            p[ ij ]      *= g->bc.gobx[i];
            v_x[ ij ]    *= g->bc.gobx[i];
            v_z[ ij ]    *= g->bc.gobx[i];
            q_x[ ij ]    *= g->bc.gobx[i];
            q_z[ ij ]    *= g->bc.gobx[i];
        }
    }
    // rightmost band, from right to left
    for ( size_t i=0; i<g->bc.np; ++i ) {
        for ( size_t j=i+1; j<g->nz2-i-1; ++j ) {
//            fprintf(fid,"%zd  %zd  %lf\n", g->nx2-i-1, j, g->bc.gobx[i]);
            
            size_t ij = (g->nx2-i-1)*g->nz2+j;
//            tau_xx[ ij ] *= g->bc.gobx[i];
//            tau_xx[ ij ] *= g->bc.gobx[i];
//            tau_xz[ ij ] *= g->bc.gobx[i];
//            p[ ij ]      *= g->bc.gobx[i];
            v_x[ ij ]    *= g->bc.gobx[i];
            v_z[ ij ]    *= g->bc.gobx[i];
            q_x[ ij ]    *= g->bc.gobx[i];
            q_z[ ij ]    *= g->bc.gobx[i];
        }
    }
//	fclose(fid);
//	exit(0);
}