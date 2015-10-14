/*
 *  pml.c
 *
 *  Created by Bernard Giroux on 11-02-11.
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

#include "pml.h"

void pml_v(double *vx_x, double *vx_z,
		   double *vz_x, double *vz_z,
		   const double *txx_x, const double *txx_z,
		   const double *tzz_x, const double *tzz_z,
		   const double *txz_x, const double *txz_z,
		   const double d[], const double dh[], const double dH[],
		   const double *txx, const double *tzz, const double *txz,
		   const double *b, const double *bm,
		   const size_t nx, const size_t nz, const size_t Npml,
		   const double dt, const double h)
{
	register double c1, c2, txx1, txx2, txz1, txz2, tzz1, tzz2;
    
	// Vx - upper left corner
	//
	size_t off1 = 0;
    
    for (size_t j=0; j<Npml; ++j) {
        txx1 = txx_x[off1 + j] + txx_z[off1 + j];
        txx2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        vx_x[ off1 + j ] *= c2;
        vx_x[ off1 + j ] += b[0] * (txx1 - txx2) / h;
        vx_x[ off1 + j ] *= c1;
    }
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txx1 = txx_x[off1 + i*Npml+j]     + txx_z[off1 + i*Npml+j];
			txx2 = txx_x[off1 + (i-1)*Npml+j] + txx_z[off1 + (i-1)*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*d[i] );
			c2 = ( 2. - dt*d[i] ) / (2.*dt);
			vx_x[ off1 + i*Npml+j ] *= c2;
			vx_x[ off1 + i*Npml+j ] += b[0] * (txx1 - txx2) / h;
			vx_x[ off1 + i*Npml+j ] *= c1;
		}
	}
	for (size_t i=0; i<Npml; ++i) {
        txz1 = txz_x[off1 + i*Npml] + txz_z[off1 + i*Npml];
		//        txz2 = -txz1;
        txz2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        vx_z[ off1 + i*Npml ] *= c2;
        vx_z[ off1 + i*Npml ] += b[0] * (txz1 - txz2) / h;
        vx_z[ off1 + i*Npml ] *= c1;
		for (size_t j=1; j<Npml; ++j) {
			txz1 = txz_x[off1 + i*Npml+j]   + txz_z[off1 + i*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j-1] + txz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[j] );
			c2 = ( 2. - dt*d[j] ) / (2.*dt);
			vx_z[ off1 + i*Npml+j ] *= c2;
			vx_z[ off1 + i*Npml+j ] += b[0] * (txz1 - txz2) / h;
			vx_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	
	// Vx - left strip
	//
	off1 += Npml*Npml;
    for (size_t j=0; j<nz; ++j) {
        txx1 = txx_x[off1 + j] + txx_z[off1 + j];
        txx2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        vx_x[ off1 + j ] *= c2;
        vx_x[ off1 + j ] += b[j] * (txx1 - txx2) / h;
        vx_x[ off1 + j ] *= c1;
    }
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<nz; ++j) {
			txx1 = txx_x[off1 + i*nz+j]     + txx_z[off1 + i*nz+j];
			txx2 = txx_x[off1 + (i-1)*nz+j] + txx_z[off1 + (i-1)*nz+j];
			c1 = (2.*dt) / ( 2. + dt*d[i] );
			c2 = ( 2. - dt*d[i] ) / (2.*dt);
			vx_x[ off1 + i*nz+j ] *= c2;
			vx_x[ off1 + i*nz+j ] += b[j] * (txx1 - txx2) / h;
			vx_x[ off1 + i*nz+j ] *= c1;
		}
	}
	size_t off2 = 0;
	for (size_t i=0; i<Npml; ++i) {
		txz1 = txz_x[off1 + i*nz]         + txz_z[off1 + i*nz];
		txz2 = txz_x[off2 + i*Npml+Npml-1] + txz_z[off2 + i*Npml+Npml-1];
		vx_z[ off1 + i*nz ] += dt * b[0] * (txz1 - txz2) / h;
		for (size_t j=1; j<nz; ++j) {
			txz1 = txz_x[off1 + i*nz+j]   + txz_z[off1 + i*nz+j];
			txz2 = txz_x[off1 + i*nz+j-1] + txz_z[off1 + i*nz+j-1];
			vx_z[ off1 + i*nz+j ] += dt * b[j] * (txz1 - txz2) / h;
		}
	}
	
	// Vx - lower left corner
	//
	off1 += Npml*nz;
    for (size_t j=0; j<Npml; ++j) {
        txx1 = txx_x[off1 + j] + txx_z[off1 + j];
        txx2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        vx_x[ off1 + j ] *= c2;
        vx_x[ off1 + j ] += b[nz-1] * (txx1 - txx2) / h;
        vx_x[ off1 + j ] *= c1;
    }
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txx1 = txx_x[off1 + i*Npml+j]     + txx_z[off1 + i*Npml+j];
			txx2 = txx_x[off1 + (i-1)*Npml+j] + txx_z[off1 + (i-1)*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*d[i] );
			c2 = ( 2. - dt*d[i] ) / (2.*dt);
			vx_x[ off1 + i*Npml+j ] *= c2;
			vx_x[ off1 + i*Npml+j ] += b[nz-1] * (txx1 - txx2) / h;
			vx_x[ off1 + i*Npml+j ] *= c1;
		}
	}
	off2 = Npml*Npml;
	for (size_t i=0; i<Npml; ++i) {
		txz1 = txz_x[off1 + i*Npml]     + txz_z[off1 + i*Npml];
		txz2 = txz_x[off2 + i*nz+nz-1] + txz_z[off2 + i*nz+nz-1];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		vx_z[ off1 + i*Npml ] *= c2;
		vx_z[ off1 + i*Npml ] += b[nz-1] * (txz1 - txz2) / h;
		vx_z[ off1 + i*Npml ] *= c1;
		for (size_t j=1; j<Npml; ++j) {
			txz1 = txz_x[off1 + i*Npml+j]   + txz_z[off1 + i*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j-1] + txz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-j] );
			c2 = ( 2. - dt*d[Npml-1-j] ) / (2.*dt);
			vx_z[ off1 + i*Npml+j ] *= c2;
			vx_z[ off1 + i*Npml+j ] += b[nz-1] * (txz1 - txz2) / h;
			vx_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	
	// Vx - upper strip
	//
	off2 = 0;
	off1 += Npml*Npml;
	for (size_t j=0; j<Npml; ++j) {
		txx1 = txx_x[off1 + j]               + txx_z[off1 + j];
		txx2 = txx_x[off2 + (Npml-1)*Npml+j] + txx_z[off2 + (Npml-1)*Npml+j];
		vx_x[ off1 + j ] += dt * b[0] * (txx1 - txx2) / h;
	}
	for (size_t i=1; i<nx; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txx1 = txx_x[off1 + i*Npml+j]     + txx_z[off1 + i*Npml+j];
			txx2 = txx_x[off1 + (i-1)*Npml+j] + txx_z[off1 + (i-1)*Npml+j];
			vx_x[ off1 + i*Npml+j ] += dt * b[i*nz] * (txx1 - txx2) / h;
		}
	}
	for (size_t i=0; i<nx; ++i) {
        txz1 = txz_x[off1 + i*Npml] + txz_z[off1 + i*Npml];
		//        txz2 = -txz1;
        txz2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        vx_z[ off1 + i*Npml ] *= c2;
        vx_z[ off1 + i*Npml ] += b[i*nz] * (txz1 - txz2) / h;
        vx_z[ off1 + i*Npml ] *= c1;
		for (size_t j=1; j<Npml; ++j) {
			txz1 = txz_x[off1 + i*Npml+j]   + txz_z[off1 + i*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j-1] + txz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[j] );
			c2 = ( 2. - dt*d[j] ) / (2.*dt);
			vx_z[ off1 + i*Npml+j ] *= c2;
			vx_z[ off1 + i*Npml+j ] += b[i*nz] * (txz1 - txz2) / h;
			vx_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	
	// Vx - lower strip
	//
	off2 = Npml*Npml + Npml*nz;
	off1 += Npml*nx;
	for (size_t j=0; j<Npml; ++j) {
		txx1 = txx_x[off1 + j]               + txx_z[off1 + j];
		txx2 = txx_x[off2 + (Npml-1)*Npml+j] + txx_z[off2 + (Npml-1)*Npml+j];
		vx_x[ off1 + j ] += dt * b[nz-1] * (txx1 - txx2) / h;
	}
	for (size_t i=1; i<nx; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txx1 = txx_x[off1 + i*Npml+j]     + txx_z[off1 + i*Npml+j];
			txx2 = txx_x[off1 + (i-1)*Npml+j] + txx_z[off1 + (i-1)*Npml+j];
			vx_x[ off1 + i*Npml+j ] += dt * b[i*nz+nz-1] * (txx1 - txx2) / h;
		}
	}
	for (size_t i=0; i<nx; ++i) {
		txz1 = txz_x[off1 + i*Npml] + txz_z[off1 + i*Npml];
		txz2 = txz[i*nz+nz-1];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		vx_z[ off1 + i*Npml ] *= c2;
		vx_z[ off1 + i*Npml ] += b[i*nz+nz-1] * (txz1 - txz2) / h;
		vx_z[ off1 + i*Npml ] *= c1;            
		for (size_t j=1; j<Npml; ++j) {
			txz1 = txz_x[off1 + i*Npml+j]   + txz_z[off1 + i*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j-1] + txz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-j] );
			c2 = ( 2. - dt*d[Npml-1-j] ) / (2.*dt);
			vx_z[ off1 + i*Npml+j ] *= c2;
			vx_z[ off1 + i*Npml+j ] += b[i*nz+nz-1] * (txz1 - txz2) / h;
			vx_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	
	// Vx - upper right corner
	//
	off2 = 2*Npml*Npml + Npml*nz;
	off1 += Npml*nx;
	for (size_t j=0; j<Npml; ++j) {  // i=0
		txx1 = txx_x[off1 + j]             + txx_z[off1 + j];
		txx2 = txx_x[off2 + (nx-1)*Npml+j] + txx_z[off2 + (nx-1)*Npml+j];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		vx_x[ off1 + j ] *= c2;
		vx_x[ off1 + j ] += b[(nx-1)*nz] * (txx1 - txx2) / h;
		vx_x[ off1 + j ] *= c1;
	}
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txx1 = txx_x[off1 + i*Npml+j]     + txx_z[off1 + i*Npml+j];
			txx2 = txx_x[off1 + (i-1)*Npml+j] + txx_z[off1 + (i-1)*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-i] );
			c2 = ( 2. - dt*d[Npml-1-i] ) / (2.*dt);
			vx_x[ off1 + i*Npml+j ] *= c2;
			vx_x[ off1 + i*Npml+j ] += b[(nx-1)*nz] * (txx1 - txx2) / h;
			vx_x[ off1 + i*Npml+j ] *= c1;
		}
	}
	for (size_t i=0; i<Npml; ++i) {
        txz1 = txz_x[off1 + i*Npml] + txz_z[off1 + i*Npml];
		//        txz2 = -txz1;
        txz2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        vx_z[ off1 + i*Npml ] *= c2;
        vx_z[ off1 + i*Npml ] += b[(nx-1)*nz] * (txz1 - txz2) / h;
        vx_z[ off1 + i*Npml ] *= c1;        
		for (size_t j=1; j<Npml; ++j) {
			txz1 = txz_x[off1 + i*Npml+j]   + txz_z[off1 + i*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j-1] + txz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[j] );
			c2 = ( 2. - dt*d[j] ) / (2.*dt);
			vx_z[ off1 + i*Npml+j ] *= c2;
			vx_z[ off1 + i*Npml+j ] += b[(nx-1)*nz] * (txz1 - txz2) / h;
			vx_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	// Vx - right strip
	//
	off1 += Npml*Npml;
	for (size_t j=0; j<nz; ++j) {  // i=0
		txx1 = txx_x[off1 + j] + txx_z[off1 + j];
		txx2 = txx[(nx-1)*nz+j];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		vx_x[ off1 + j ] *= c2;
		vx_x[ off1 + j ] += b[(nx-1)*nz+j] * (txx1 - txx2) / h;
		vx_x[ off1 + j ] *= c1;
	}
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<nz; ++j) {
			txx1 = txx_x[off1 + i*nz+j]     + txx_z[off1 + i*nz+j];
			txx2 = txx_x[off1 + (i-1)*nz+j] + txx_z[off1 + (i-1)*nz+j];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-i] );
			c2 = ( 2. - dt*d[Npml-1-i] ) / (2.*dt);
			vx_x[ off1 + i*nz+j ] *= c2;
			vx_x[ off1 + i*nz+j ] += b[(nx-1)*nz+j] * (txx1 - txx2) / h;
			vx_x[ off1 + i*nz+j ] *= c1;
		}
	}
	off2 = 2*Npml*Npml + Npml*nz + 2*Npml*nx;
	for (size_t i=0; i<Npml; ++i) {
		txz1 = txz_x[off1 + i*nz]          + txz_z[off1 + i*nz];
		txz2 = txz_x[off2 + i*Npml+Npml-1] + txz_z[off2 + i*Npml+Npml-1];
		vx_z[ off1 + i*nz ] += dt * b[(nx-1)*nz] * (txz1 - txz2) / h;
		for (size_t j=1; j<nz; ++j) {
			txz1 = txz_x[off1 + i*nz+j]   + txz_z[off1 + i*nz+j];
			txz2 = txz_x[off1 + i*nz+j-1] + txz_z[off1 + i*nz+j-1];
			vx_z[ off1 + i*nz+j ] += dt * b[(nx-1)*nz+j] * (txz1 - txz2) / h;
		}
	}
	// Vx - lower right corner
	//
	off2 = 2*Npml*Npml + Npml*nz + Npml*nx;
	off1 += Npml*nz;
	for (size_t j=0; j<Npml; ++j) {  // i=0
		txx1 = txx_x[off1 + j]             + txx_z[off1 + j];
		txx2 = txx_x[off2 + (nx-1)*Npml+j] + txx_z[off2 + (nx-1)*Npml+j];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		vx_x[ off1 + j ] *= c2;
		vx_x[ off1 + j ] += b[nx*nz-1] * (txx1 - txx2) / h;
		vx_x[ off1 + j ] *= c1;
	}
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txx1 = txx_x[off1 + i*Npml+j]     + txx_z[off1 + i*Npml+j];
			txx2 = txx_x[off1 + (i-1)*Npml+j] + txx_z[off1 + (i-1)*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-i] );
			c2 = ( 2. - dt*d[Npml-1-i] ) / (2.*dt);
			vx_x[ off1 + i*Npml+j ] *= c2;
			vx_x[ off1 + i*Npml+j ] += b[nx*nz-1] * (txx1 - txx2) / h;
			vx_x[ off1 + i*Npml+j ] *= c1;
		}
	}
	off2 = 3*Npml*Npml + Npml*nz + 2*Npml*nx;
	for (size_t i=0; i<Npml; ++i) {
		txz1 = txz_x[off1 + i*Npml]    + txz_z[off1 + i*Npml];
		txz2 = txz_x[off2 + i*nz+nz-1] + txz_z[off2 + i*nz+nz-1];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		vx_z[ off1 + i*Npml ] *= c2;
		vx_z[ off1 + i*Npml ] += b[nx*nz-1] * (txz1 - txz2) / h;
		vx_z[ off1 + i*Npml ] *= c1;
		for (size_t j=1; j<Npml; ++j) {
			txz1 = txz_x[off1 + i*Npml+j]   + txz_z[off1 + i*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j-1] + txz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-j] );
			c2 = ( 2. - dt*d[Npml-1-j] ) / (2.*dt);
			vx_z[ off1 + i*Npml+j ] *= c2;
			vx_z[ off1 + i*Npml+j ] += b[nx*nz-1] * (txz1 - txz2) / h;
			vx_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	
	
	
	
	// Vz - upper left corner
	//
	off1 = 0;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txz1 = txz_x[off1 + (i+1)*Npml+j] + txz_z[off1 + (i+1)*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j]     + txz_z[off1 + i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dh[i] );
			c2 = ( 2. - dt*dh[i] ) / (2.*dt);
			vz_x[ off1 + i*Npml+j ] *= c2;
			vz_x[ off1 + i*Npml+j ] += bm[0] * (txz1 - txz2) / h;
			vz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
	off2 = 2*Npml*Npml + Npml*nz;
	for (size_t j=0; j<Npml; ++j) {
		txz1 = txz_x[off2 + j]                 + txz_z[off2 + j];
		txz2 = txz_x[off1 + (Npml-1)*Npml+j]   + txz_z[off1 + (Npml-1)*Npml+j];
		c1 = (2.*dt) / ( 2. + dt*dh[Npml-1] );
		c2 = ( 2. - dt*dh[Npml-1] ) / (2.*dt);
		vz_x[ off1 + (Npml-1)*Npml+j ] *= c2;
		vz_x[ off1 + (Npml-1)*Npml+j ] += bm[0] * (txz1 - txz2) / h;
		vz_x[ off1 + (Npml-1)*Npml+j ] *= c1;
	}
	off2 = Npml*Npml;
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			tzz1 = tzz_x[ off1 + i*Npml+j+1 ] + tzz_z[ off1 + i*Npml+j+1 ];
			tzz2 = tzz_x[ off1 + i*Npml+j ]   + tzz_z[ off1 + i*Npml+j ];
			c1 = (2.*dt) / ( 2. + dt*dh[j] );
			c2 = ( 2. - dt*dh[j] ) / (2.*dt);
			vz_z[ off1 + i*Npml+j ] *= c2;
			vz_z[ off1 + i*Npml+j ] += bm[0] * (tzz1 - tzz2) / h;
			vz_z[ off1 + i*Npml+j ] *= c1;
		}
		tzz1 = tzz_x[ off2 + i*nz ]          + tzz_z[ off2 + i*nz ];
		tzz2 = tzz_x[ off1 + i*Npml+Npml-1 ] + tzz_z[ off1 + i*Npml+Npml-1 ];
		c1 = (2.*dt) / ( 2. + dt*dh[Npml-1] );
		c2 = ( 2. - dt*dh[Npml-1] ) / (2.*dt);
		vz_z[ off1 + i*Npml+Npml-1 ] *= c2;
		vz_z[ off1 + i*Npml+Npml-1 ] += bm[0] * (tzz1 - tzz2) / h;
		vz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}
	// Vz - left strip
	//
	off1 += Npml*Npml;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<nz; ++j) {
			txz1 = txz_x[off1 + (i+1)*nz+j] + txz_z[off1 + (i+1)*nz+j];
			txz2 = txz_x[off1 + i*nz+j]     + txz_z[off1 + i*nz+j];
			c1 = (2.*dt) / ( 2. + dt*dh[i] );
			c2 = ( 2. - dt*dh[i] ) / (2.*dt);
			vz_x[ off1 + i*nz+j ] *= c2;
			vz_x[ off1 + i*nz+j ] += bm[j] * (txz1 - txz2) / h;
			vz_x[ off1 + i*nz+j ] *= c1;
		}
	}
	for (size_t j=0; j<nz; ++j) {
		txz1 = txz[j];
		txz2 = txz_x[off1 + (Npml-1)*nz+j] + txz_z[off1 + (Npml-1)*nz+j];
		c1 = (2.*dt) / ( 2. + dt*dh[(Npml-1)] );
		c2 = ( 2. - dt*dh[(Npml-1)] ) / (2.*dt);
		vz_x[ off1 + (Npml-1)*nz+j ] *= c2;
		vz_x[ off1 + (Npml-1)*nz+j ] += bm[j] * (txz1 - txz2) / h;
		vz_x[ off1 + (Npml-1)*nz+j ] *= c1;
	}
	off2 = Npml*Npml + Npml*nz;
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<nz-1; ++j) {
			tzz1 = tzz_x[ off1 + i*nz+j+1 ] + tzz_z[ off1 + i*nz+j+1 ];
			tzz2 = tzz_x[ off1 + i*nz+j ]   + tzz_z[ off1 + i*nz+j ];
			vz_z[ off1 + i*nz+j ] += dt * bm[j] * (tzz1 - tzz2) / h;
		}
		tzz1 = tzz_x[ off2 + i*Npml ]    + tzz_z[ off2 + i*Npml ];
		tzz2 = tzz_x[ off1 + i*nz+nz-1 ] + tzz_z[ off1 + i*nz+nz-1 ];
		vz_z[ off1 + i*nz+nz-1 ] += dt * bm[nz-1] * (tzz1 - tzz2) / h;
	}                
	// Vz - lower left corner
	//
	off1 += Npml*nz;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txz1 = txz_x[off1 + (i+1)*Npml+j] + txz_z[off1 + (i+1)*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j]     + txz_z[off1 + i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dh[i] );
			c2 = ( 2. - dt*dh[i] ) / (2.*dt);
			vz_x[ off1 + i*Npml+j ] *= c2;
			vz_x[ off1 + i*Npml+j ] += bm[nz-1] * (txz1 - txz2) / h;
			vz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
	off2 = 2*Npml*Npml + Npml*nz + Npml*nx;
	for (size_t j=0; j<Npml; ++j) {
		txz1 = txz_x[off2 + j]                 + txz_z[off2 + j];
		txz2 = txz_x[off1 + (Npml-1)*Npml+j]   + txz_z[off1 + (Npml-1)*Npml+j];
		c1 = (2.*dt) / ( 2. + dt*dh[Npml-1] );
		c2 = ( 2. - dt*dh[Npml-1] ) / (2.*dt);
		vz_x[ off1 + (Npml-1)*Npml+j ] *= c2;
		vz_x[ off1 + (Npml-1)*Npml+j ] += bm[nz-1] * (txz1 - txz2) / h;
		vz_x[ off1 + (Npml-1)*Npml+j ] *= c1;
	}
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			tzz1 = tzz_x[ off1 + i*Npml+j+1 ] + tzz_z[ off1 + i*Npml+j+1 ];
			tzz2 = tzz_x[ off1 + i*Npml+j ]   + tzz_z[ off1 + i*Npml+j ];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-j] );
			c2 = ( 2. - dt*dH[Npml-1-j] ) / (2.*dt);
			vz_z[ off1 + i*Npml+j ] *= c2;
			vz_z[ off1 + i*Npml+j ] += bm[nz-1] * (tzz1 - tzz2) / h;
			vz_z[ off1 + i*Npml+j ] *= c1;
		}
        tzz1 = 0.;
        tzz2 = tzz_x[ off1 + i*Npml+Npml-1 ]   + tzz_z[ off1 + i*Npml+Npml-1 ];
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        vz_z[ off1 + i*Npml+Npml-1 ] *= c2;
        vz_z[ off1 + i*Npml+Npml-1 ] += bm[nz-1] * (tzz1 - tzz2) / h;
        vz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}  
	
	// Vz - upper strip
	//
	off1 += Npml*Npml;
	for (size_t i=0; i<nx-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txz1 = txz_x[off1 + (i+1)*Npml+j] + txz_z[off1 + (i+1)*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j]     + txz_z[off1 + i*Npml+j];
			vz_x[ off1 + i*Npml+j ] += dt * bm[i*nz] * (txz1 - txz2) / h;
		}
	}
	off2 = 2*Npml*Npml + Npml*nz + 2*Npml*nx;
	for (size_t j=0; j<Npml; ++j) {
		txz1 = txz_x[off2 + j]               + txz_z[off2 + j];
		txz2 = txz_x[off1 + (nx-1)*Npml+j]   + txz_z[off1 + (nx-1)*Npml+j];
		vz_x[ off1 + (nx-1)*Npml+j ] += dt * bm[(nx-1)*nz] * (txz1 - txz2) / h;
	}
	for (size_t i=0; i<nx; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			tzz1 = tzz_x[ off1 + i*Npml+j+1 ] + tzz_z[ off1 + i*Npml+j+1 ];
			tzz2 = tzz_x[ off1 + i*Npml+j ]   + tzz_z[ off1 + i*Npml+j ];
			c1 = (2.*dt) / ( 2. + dt*dh[j] );
			c2 = ( 2. - dt*dh[j] ) / (2.*dt);
			vz_z[ off1 + i*Npml+j ] *= c2;
			vz_z[ off1 + i*Npml+j ] += bm[i*nz] * (tzz1 - tzz2) / h;
			vz_z[ off1 + i*Npml+j ] *= c1;
		}
		tzz1 = tzz[ i*nz ];
		tzz2 = tzz_x[ off1 + i*Npml+Npml-1 ] + tzz_z[ off1 + i*Npml+Npml-1 ];
		c1 = (2.*dt) / ( 2. + dt*dh[Npml-1] );
		c2 = ( 2. - dt*dh[Npml-1] ) / (2.*dt);
		vz_z[ off1 + i*Npml+Npml-1 ] *= c2;
		vz_z[ off1 + i*Npml+Npml-1 ] += bm[i*nz] * (tzz1 - tzz2) / h;
		vz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}
	
	// Vz - lower strip
	//
	off1 += Npml*nx;
	for (size_t i=0; i<nx-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txz1 = txz_x[off1 + (i+1)*Npml+j] + txz_z[off1 + (i+1)*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j]     + txz_z[off1 + i*Npml+j];
			vz_x[ off1 + i*Npml+j ] += dt * bm[i*nz+nz-1] * (txz1 - txz2) / h;
		}
	}
	off2 = 3*Npml*Npml + 2*Npml*nz + 2*Npml*nx;
	for (size_t j=0; j<Npml; ++j) {
		txz1 = txz_x[off2 + j]             + txz_z[off2 + j];
		txz2 = txz_x[off1 + (nx-1)*Npml+j] + txz_z[off1 + (nx-1)*Npml+j];
		vz_x[ off1 + (nx-1)*Npml+j ] += dt * bm[(nx-1)*nz+nz-1] * (txz1 - txz2) / h;
	}
	for (size_t i=0; i<nx; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			tzz1 = tzz_x[ off1 + i*Npml+j+1 ] + tzz_z[ off1 + i*Npml+j+1 ];
			tzz2 = tzz_x[ off1 + i*Npml+j ]   + tzz_z[ off1 + i*Npml+j ];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-j] );
			c2 = ( 2. - dt*dH[Npml-1-j] ) / (2.*dt);
			vz_z[ off1 + i*Npml+j ] *= c2;
			vz_z[ off1 + i*Npml+j ] += bm[i*nz+nz-1] * (tzz1 - tzz2) / h;
			vz_z[ off1 + i*Npml+j ] *= c1;
		}
        tzz1 = 0.;
        tzz2 = tzz_x[ off1 + i*Npml+Npml-1 ]   + tzz_z[ off1 + i*Npml+Npml-1 ];
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        vz_z[ off1 + i*Npml+Npml-1 ] *= c2;
        vz_z[ off1 + i*Npml+Npml-1 ] += bm[i*nz+nz-1] * (tzz1 - tzz2) / h;
        vz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}
	// Vz - upper right corner
	//
	off1 += Npml*nx;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txz1 = txz_x[off1 + (i+1)*Npml+j] + txz_z[off1 + (i+1)*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j]     + txz_z[off1 + i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-i] );
			c2 = ( 2. - dt*dH[Npml-1-i] ) / (2.*dt);
			vz_x[ off1 + i*Npml+j ] *= c2;
			vz_x[ off1 + i*Npml+j ] += bm[(nx-1)*nz] * (txz1 - txz2) / h;
			vz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
    for (size_t j=0; j<Npml; ++j) {
        txz1 = 0.;
        txz2 = txz_x[off1 + (Npml-1)*Npml+j] + txz_z[off1 + (Npml-1)*Npml+j];
		//		txz1 = -txz2;
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        vz_x[ off1 + (Npml-1)*Npml+j ] *= c2;
        vz_x[ off1 + (Npml-1)*Npml+j ] += bm[(nx-1)*nz] * (txz1 - txz2) / h;
        vz_x[ off1 + (Npml-1)*Npml+j ] *= c1;
    }
	off2 = 3*Npml*Npml + Npml*nz + 2*Npml*nx;
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			tzz1 = tzz_x[ off1 + i*Npml+j+1 ] + tzz_z[ off1 + i*Npml+j+1 ];
			tzz2 = tzz_x[ off1 + i*Npml+j ]   + tzz_z[ off1 + i*Npml+j ];
			c1 = (2.*dt) / ( 2. + dt*dh[j] );
			c2 = ( 2. - dt*dh[j] ) / (2.*dt);
			vz_z[ off1 + i*Npml+j ] *= c2;
			vz_z[ off1 + i*Npml+j ] += bm[(nx-1)*nz] * (tzz1 - tzz2) / h;
			vz_z[ off1 + i*Npml+j ] *= c1;
		}
		tzz1 = tzz_x[ off2 + i*nz ]          + tzz_z[ off2 + i*nz ];
		tzz2 = tzz_x[ off1 + i*Npml+Npml-1 ] + tzz_z[ off1 + i*Npml+Npml-1 ];
		c1 = (2.*dt) / ( 2. + dt*dh[Npml-1] );
		c2 = ( 2. - dt*dh[Npml-1] ) / (2.*dt);
		vz_z[ off1 + i*Npml+Npml-1 ] *= c2;
		vz_z[ off1 + i*Npml+Npml-1 ] += bm[(nx-1)*nz] * (tzz1 - tzz2) / h;
		vz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}
	
	// Vz - right strip
	//
	off1 += Npml*Npml;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<nz; ++j) {
			txz1 = txz_x[off1 + (i+1)*nz+j] + txz_z[off1 + (i+1)*nz+j];
			txz2 = txz_x[off1 + i*nz+j]     + txz_z[off1 + i*nz+j];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-i] );
			c2 = ( 2. - dt*dH[Npml-1-i] ) / (2.*dt);
			vz_x[ off1 + i*nz+j ] *= c2;
			vz_x[ off1 + i*nz+j ] += bm[(nx-1)*nz+j] * (txz1 - txz2) / h;
			vz_x[ off1 + i*nz+j ] *= c1;
		}
	}
    for (size_t j=0; j<nz; ++j) {
        txz1 = 0.;
        txz2 = txz_x[off1 + (Npml-1)*nz+j] + txz_z[off1 + (Npml-1)*nz+j];
		//		txz1 = -txz2;
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        vz_x[ off1 + (Npml-1)*nz+j ] *= c2;
        vz_x[ off1 + (Npml-1)*nz+j ] += bm[(nx-1)*nz+j] * (txz1 - txz2) / h;
        vz_x[ off1 + (Npml-1)*nz+j ] *= c1;
    }
	off2 = 3*Npml*Npml + 2*Npml*nz + 2*Npml*nx;
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<nz-1; ++j) {
			tzz1 = tzz_x[ off1 + i*nz+j+1 ] + tzz_z[ off1 + i*nz+j+1 ];
			tzz2 = tzz_x[ off1 + i*nz+j ]   + tzz_z[ off1 + i*nz+j ];
			vz_z[ off1 + i*nz+j ] += dt * bm[(nx-1)*nz+j] * (tzz1 - tzz2) / h;
		}
		tzz1 = tzz_x[ off2 + i*Npml ]    + tzz_z[ off2 + i*Npml ];
		tzz2 = tzz_x[ off1 + i*nz+nz-1 ] + tzz_z[ off1 + i*nz+nz-1 ];
		vz_z[ off1 + i*nz+nz-1 ] += dt * bm[nx*nz-1] * (tzz1 - tzz2) / h;
	}
	
	// Vz - lower right corner
	//
	off1 += Npml*nz;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			txz1 = txz_x[off1 + (i+1)*Npml+j] + txz_z[off1 + (i+1)*Npml+j];
			txz2 = txz_x[off1 + i*Npml+j]     + txz_z[off1 + i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-i] );
			c2 = ( 2. - dt*dH[Npml-1-i] ) / (2.*dt);
			vz_x[ off1 + i*Npml+j ] *= c2;
			vz_x[ off1 + i*Npml+j ] += bm[nx*nz-1] * (txz1 - txz2) / h;
			vz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
    for (size_t j=0; j<Npml; ++j) {
        txz1 = 0.;
        txz2 = txz_x[off1 + (Npml-1)*Npml+j] + txz_z[off1 + (Npml-1)*Npml+j];
		//		txz1 = -txz2;
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        vz_x[ off1 + (Npml-1)*Npml+j ] *= c2;
        vz_x[ off1 + (Npml-1)*Npml+j ] += bm[nx*nz-1] * (txz1 - txz2) / h;
        vz_x[ off1 + (Npml-1)*Npml+j ] *= c1;
    }
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			tzz1 = tzz_x[ off1 + i*Npml+j+1 ] + tzz_z[ off1 + i*Npml+j+1 ];
			tzz2 = tzz_x[ off1 + i*Npml+j ]   + tzz_z[ off1 + i*Npml+j ];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-j] );
			c2 = ( 2. - dt*dH[Npml-1-j] ) / (2.*dt);
			vz_z[ off1 + i*Npml+j ] *= c2;
			vz_z[ off1 + i*Npml+j ] += bm[nx*nz-1] * (tzz1 - tzz2) / h;
			vz_z[ off1 + i*Npml+j ] *= c1;
		}
        tzz1 = 0.;
        tzz2 = tzz_x[ off1 + i*Npml+Npml-1 ] + tzz_z[ off1 + i*Npml+Npml-1 ];
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        vz_z[ off1 + i*Npml+Npml-1 ] *= c2;
        vz_z[ off1 + i*Npml+Npml-1 ] += bm[nx*nz-1] * (tzz1 - tzz2) / h;
        vz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}
}



void pml_t(const double *vx_x, const double *vx_z,
		   const double *vz_x, const double *vz_z,
		   double *txx_x, double *txx_z,
		   double *tzz_x, double *tzz_z,
		   double *txz_x, double *txz_z,
		   const double d[], const double dh[], const double dH[],
		   const double *vx, const double *vz,
		   const double *la, const double *l2m, const double *mu,
		   const size_t nx, const size_t nz, const size_t Npml,
		   const double dt, const double h)
{
    
    register double c1, c2, vx1, vx2, vz1, vz2;
    
	// txx & tzz - upper left corner
	//
	size_t off1 = 0;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vx1 = vx_x[off1 + (i+1)*Npml+j] + vx_z[off1 + (i+1)*Npml+j];
			vx2 = vx_x[off1 +     i*Npml+j] + vx_z[off1 +     i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dh[i] );
			c2 = ( 2. - dt*dh[i] ) / (2.*dt);
			txx_x[ off1 + i*Npml+j ] *= c2;
			txx_x[ off1 + i*Npml+j ] += l2m[0] * (vx1 - vx2) / h;
			txx_x[ off1 + i*Npml+j ] *= c1;
			tzz_x[ off1 + i*Npml+j ] *= c2;
			tzz_x[ off1 + i*Npml+j ] += la[0] * (vx1 - vx2) / h;
			tzz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
    size_t off2 = 2*Npml*Npml + Npml*nz;
	for (size_t j=0; j<Npml; ++j) {
		vx1 = vx_x[off2 + j]               + vx_z[off2 + j];
		vx2 = vx_x[off1 + (Npml-1)*Npml+j] + vx_z[off1 + (Npml-1)*Npml+j];
		c1 = (2.*dt) / ( 2. + dt*dh[(Npml-1)] );
		c2 = ( 2. - dt*dh[(Npml-1)] ) / (2.*dt);
		txx_x[ off1 + (Npml-1)*Npml+j ] *= c2;
		txx_x[ off1 + (Npml-1)*Npml+j ] += l2m[0] * (vx1 - vx2) / h;
		txx_x[ off1 + (Npml-1)*Npml+j ] *= c1;
		tzz_x[ off1 + (Npml-1)*Npml+j ] *= c2;
		tzz_x[ off1 + (Npml-1)*Npml+j ] += la[0] * (vx1 - vx2) / h;
		tzz_x[ off1 + (Npml-1)*Npml+j ] *= c1;
	}
	for (size_t i=0; i<Npml; ++i) {
        vz1 = vz_x[off1 + i*Npml] + vz_z[off1 + i*Npml];
        vz2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        txx_z[ off1 + i*Npml ] *= c2;
        txx_z[ off1 + i*Npml ] += la[0] * (vz1 - vz2) / h;
        txx_z[ off1 + i*Npml ] *= c1;
        tzz_z[ off1 + i*Npml ] *= c2;
        tzz_z[ off1 + i*Npml ] += l2m[0] * (vz1 - vz2) / h;
        tzz_z[ off1 + i*Npml ] *= c1;
		for (size_t j=1; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]   + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + i*Npml+j-1] + vz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[j] );
			c2 = ( 2. - dt*d[j] ) / (2.*dt);
			txx_z[ off1 + i*Npml+j ] *= c2;
			txx_z[ off1 + i*Npml+j ] += la[0] * (vz1 - vz2) / h;
			txx_z[ off1 + i*Npml+j ] *= c1;
			tzz_z[ off1 + i*Npml+j ] *= c2;
			tzz_z[ off1 + i*Npml+j ] += l2m[0] * (vz1 - vz2) / h;
			tzz_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	// txx & tzz - left strip
	//
	off1 += Npml*Npml;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<nz; ++j) {
			vx1 = vx_x[off1 + (i+1)*nz+j] + vx_z[off1 + (i+1)*nz+j];
			vx2 = vx_x[off1 +     i*nz+j] + vx_z[off1 +     i*nz+j];
			c1 = (2.*dt) / ( 2. + dt*dh[i] );
			c2 = ( 2. - dt*dh[i] ) / (2.*dt);
			txx_x[ off1 + i*nz+j ] *= c2;
			txx_x[ off1 + i*nz+j ] += l2m[j] * (vx1 - vx2) / h;
			txx_x[ off1 + i*nz+j ] *= c1;
			tzz_x[ off1 + i*nz+j ] *= c2;
			tzz_x[ off1 + i*nz+j ] += la[j] * (vx1 - vx2) / h;
			tzz_x[ off1 + i*nz+j ] *= c1;
		}
	}
	for (size_t j=0; j<nz; ++j) {
		vx1 = vx[j];
		vx2 = vx_x[off1 + (Npml-1)*nz+j] + vx_z[off1 + (Npml-1)*nz+j];
		c1 = (2.*dt) / ( 2. + dt*dh[(Npml-1)] );
		c2 = ( 2. - dt*dh[(Npml-1)] ) / (2.*dt);
		txx_x[ off1 + (Npml-1)*nz+j ] *= c2;
		txx_x[ off1 + (Npml-1)*nz+j ] += l2m[j] * (vx1 - vx2) / h;
		txx_x[ off1 + (Npml-1)*nz+j ] *= c1;
		tzz_x[ off1 + (Npml-1)*nz+j ] *= c2;
		tzz_x[ off1 + (Npml-1)*nz+j ] += la[j] * (vx1 - vx2) / h;
		tzz_x[ off1 + (Npml-1)*nz+j ] *= c1;
	}
	off2 = 0;
	for (size_t i=0; i<Npml; ++i) {
		vz1 = vz_x[off1 + i*nz]          + vz_z[off1 + i*nz];
		vz2 = vz_x[off2 + i*Npml+Npml-1] + vz_z[off2 + i*Npml+Npml-1];
		txx_z[ off1 + i*nz ] += dt * la[0] * (vz1 - vz2) / h;
		tzz_z[ off1 + i*nz ] += dt * l2m[0] * (vz1 - vz2) / h;
		for (size_t j=1; j<nz; ++j) {
			vz1 = vz_x[off1 + i*nz+j]   + vz_z[off1 + i*nz+j];
			vz2 = vz_x[off1 + i*nz+j-1] + vz_z[off1 + i*nz+j-1];
			txx_z[ off1 + i*nz+j ] += dt * la[j] * (vz1 - vz2) / h;
			tzz_z[ off1 + i*nz+j ] += dt * l2m[j] * (vz1 - vz2) / h;
		}
	}
	
	// txx & tzz - lower left corner
	//
	off1 += Npml*nz;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vx1 = vx_x[off1 + (i+1)*Npml+j] + vx_z[off1 + (i+1)*Npml+j];
			vx2 = vx_x[off1 +     i*Npml+j] + vx_z[off1 +     i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dh[i] );
			c2 = ( 2. - dt*dh[i] ) / (2.*dt);
			txx_x[ off1 + i*Npml+j ] *= c2;
			txx_x[ off1 + i*Npml+j ] += l2m[nz-1] * (vx1 - vx2) / h;
			txx_x[ off1 + i*Npml+j ] *= c1;
			tzz_x[ off1 + i*Npml+j ] *= c2;
			tzz_x[ off1 + i*Npml+j ] += la[nz-1] * (vx1 - vx2) / h;
			tzz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
    off2 = 2*Npml*Npml + Npml*nz + Npml*nx;
	for (size_t j=0; j<Npml; ++j) {
		vx1 = vx_x[off2 + j]               + vx_z[off2 + j];
		vx2 = vx_x[off1 + (Npml-1)*Npml+j] + vx_z[off1 + (Npml-1)*Npml+j];
		c1 = (2.*dt) / ( 2. + dt*dh[(Npml-1)] );
		c2 = ( 2. - dt*dh[(Npml-1)] ) / (2.*dt);
		txx_x[ off1 + (Npml-1)*Npml+j ] *= c2;
		txx_x[ off1 + (Npml-1)*Npml+j ] += l2m[nz-1] * (vx1 - vx2) / h;
		txx_x[ off1 + (Npml-1)*Npml+j ] *= c1;
		tzz_x[ off1 + (Npml-1)*Npml+j ] *= c2;
		tzz_x[ off1 + (Npml-1)*Npml+j ] += la[nz-1] * (vx1 - vx2) / h;
		tzz_x[ off1 + (Npml-1)*Npml+j ] *= c1;
	}
	off2 = Npml*Npml;
	for (size_t i=0; i<Npml; ++i) {
		vz1 = vz_x[off1 + i*Npml]    + vz_z[off1 + i*Npml];
		vz2 = vz_x[off2 + i*nz+nz-1] + vz_z[off2 + i*nz+nz-1];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		txx_z[ off1 + i*Npml ] *= c2;
		txx_z[ off1 + i*Npml ] += la[nz-1] * (vz1 - vz2) / h;
		txx_z[ off1 + i*Npml ] *= c1;
		tzz_z[ off1 + i*Npml ] *= c2;
		tzz_z[ off1 + i*Npml ] += l2m[nz-1] * (vz1 - vz2) / h;
		tzz_z[ off1 + i*Npml ] *= c1;
		for (size_t j=1; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]   + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + i*Npml+j-1] + vz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-j] );
			c2 = ( 2. - dt*d[Npml-1-j] ) / (2.*dt);
			txx_z[ off1 + i*Npml+j ] *= c2;
			txx_z[ off1 + i*Npml+j ] += la[nz-1] * (vz1 - vz2) / h;
			txx_z[ off1 + i*Npml+j ] *= c1;
			tzz_z[ off1 + i*Npml+j ] *= c2;
			tzz_z[ off1 + i*Npml+j ] += l2m[nz-1] * (vz1 - vz2) / h;
			tzz_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	
	// txx & tzz - upper strip
	//
	off1 += Npml*Npml;
	for (size_t i=0; i<nx-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vx1 = vx_x[off1 + (i+1)*Npml+j] + vx_z[off1 + (i+1)*Npml+j];
			vx2 = vx_x[off1 +     i*Npml+j] + vx_z[off1 +     i*Npml+j];
			txx_x[ off1 + i*Npml+j ] += dt * l2m[i*nz] * (vx1 - vx2) / h;
			tzz_x[ off1 + i*Npml+j ] += dt *  la[i*nz] * (vx1 - vx2) / h;
		}
	}
	off2 = 2*Npml*Npml + Npml*nz + 2*Npml*nx;
	for (size_t j=0; j<Npml; ++j) {
		vx1 = vx_x[off2 + j]             + vx_z[off2 + j];
		vx2 = vx_x[off1 + (nx-1)*Npml+j] + vx_z[off1 + (nx-1)*Npml+j];
		txx_x[ off1 + (nx-1)*Npml+j ] += dt * l2m[(nx-1)*nz] * (vx1 - vx2) / h;
		tzz_x[ off1 + (nx-1)*Npml+j ] += dt *  la[(nx-1)*nz] * (vx1 - vx2) / h;
	}
	for (size_t i=0; i<nx; ++i) {
        vz1 = vz_x[off1 + i*Npml] + vz_z[off1 + i*Npml];
        vz2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        txx_z[ off1 + i*Npml ] *= c2;
        txx_z[ off1 + i*Npml ] += la[i*nz] * (vz1 - vz2) / h;
        txx_z[ off1 + i*Npml ] *= c1;
        tzz_z[ off1 + i*Npml ] *= c2;
        tzz_z[ off1 + i*Npml ] += l2m[i*nz] * (vz1 - vz2) / h;
        tzz_z[ off1 + i*Npml ] *= c1;
		for (size_t j=1; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]   + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + i*Npml+j-1] + vz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[j] );
			c2 = ( 2. - dt*d[j] ) / (2.*dt);
			txx_z[ off1 + i*Npml+j ] *= c2;
			txx_z[ off1 + i*Npml+j ] += la[i*nz] * (vz1 - vz2) / h;
			txx_z[ off1 + i*Npml+j ] *= c1;
			tzz_z[ off1 + i*Npml+j ] *= c2;
			tzz_z[ off1 + i*Npml+j ] += l2m[i*nz] * (vz1 - vz2) / h;
			tzz_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	
	// txx & tzz - lower strip
	//
	off1 += Npml*nx;
	for (size_t i=0; i<nx-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vx1 = vx_x[off1 + (i+1)*Npml+j] + vx_z[off1 + (i+1)*Npml+j];
			vx2 = vx_x[off1 +     i*Npml+j] + vx_z[off1 +     i*Npml+j];
			txx_x[ off1 + i*Npml+j ] += dt * l2m[i*nz+nz-1] * (vx1 - vx2) / h;
			tzz_x[ off1 + i*Npml+j ] += dt *  la[i*nz+nz-1] * (vx1 - vx2) / h;
		}
	}
	off2 = 3*Npml*Npml + 2*Npml*nz + 2*Npml*nx;
	for (size_t j=0; j<Npml; ++j) {
		vx1 = vx_x[off2 + j]             + vx_z[off2 + j];
		vx2 = vx_x[off1 + (nx-1)*Npml+j] + vx_z[off1 + (nx-1)*Npml+j];
		txx_x[ off1 + (nx-1)*Npml+j ] += dt * l2m[nx*nz-1] * (vx1 - vx2) / h;
		tzz_x[ off1 + (nx-1)*Npml+j ] += dt *  la[nx*nz-1] * (vx1 - vx2) / h;
	}
	for (size_t i=0; i<nx; ++i) {
		vz1 = vz_x[off1 + i*Npml] + vz_z[off1 + i*Npml];
		vz2 = vz[i*nz+nz-1];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		txx_z[ off1 + i*Npml ] *= c2;
		txx_z[ off1 + i*Npml ] += la[i*nz+nz-1] * (vz1 - vz2) / h;
		txx_z[ off1 + i*Npml ] *= c1;
		tzz_z[ off1 + i*Npml ] *= c2;
		tzz_z[ off1 + i*Npml ] += l2m[i*nz+nz-1] * (vz1 - vz2) / h;
		tzz_z[ off1 + i*Npml ] *= c1;
		for (size_t j=1; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]   + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + i*Npml+j-1] + vz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-j] );
			c2 = ( 2. - dt*d[Npml-1-j] ) / (2.*dt);
			txx_z[ off1 + i*Npml+j ] *= c2;
			txx_z[ off1 + i*Npml+j ] += la[i*nz+nz-1] * (vz1 - vz2) / h;
			txx_z[ off1 + i*Npml+j ] *= c1;
			tzz_z[ off1 + i*Npml+j ] *= c2;
			tzz_z[ off1 + i*Npml+j ] += l2m[i*nz+nz-1] * (vz1 - vz2) / h;
			tzz_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	
	// txx & tzz - upper right corner
	//
	off1 += Npml*nx;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vx1 = vx_x[off1 + (i+1)*Npml+j] + vx_z[off1 + (i+1)*Npml+j];
			vx2 = vx_x[off1 +     i*Npml+j] + vx_z[off1 +     i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-i] );
			c2 = ( 2. - dt*dH[Npml-1-i] ) / (2.*dt);
			txx_x[ off1 + i*Npml+j ] *= c2;
			txx_x[ off1 + i*Npml+j ] += l2m[(nx-1)*nz] * (vx1 - vx2) / h;
			txx_x[ off1 + i*Npml+j ] *= c1;
			tzz_x[ off1 + i*Npml+j ] *= c2;
			tzz_x[ off1 + i*Npml+j ] += la[(nx-1)*nz] * (vx1 - vx2) / h;
			tzz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
    for (size_t j=0; j<Npml; ++j) {
        vx1 = 0.;
        vx2 = vx_x[off1 + (Npml-1)*Npml+j] + vx_z[off1 + (Npml-1)*Npml+j];
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        txx_x[ off1 + (Npml-1)*Npml+j ] *= c2;
        txx_x[ off1 + (Npml-1)*Npml+j ] += l2m[(nx-1)*nz] * (vx1 - vx2) / h;
        txx_x[ off1 + (Npml-1)*Npml+j ] *= c1;
        tzz_x[ off1 + (Npml-1)*Npml+j ] *= c2;
        tzz_x[ off1 + (Npml-1)*Npml+j ] += la[(nx-1)*nz] * (vx1 - vx2) / h;
        tzz_x[ off1 + (Npml-1)*Npml+j ] *= c1;
    }
	for (size_t i=0; i<Npml; ++i) {
        vz1 = vz_x[off1 + i*Npml] + vz_z[off1 + i*Npml];
        vz2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        txx_z[ off1 + i*Npml ] *= c2;
        txx_z[ off1 + i*Npml ] += la[(nx-1)*nz] * (vz1 - vz2) / h;
        txx_z[ off1 + i*Npml ] *= c1;
        tzz_z[ off1 + i*Npml ] *= c2;
        tzz_z[ off1 + i*Npml ] += l2m[(nx-1)*nz] * (vz1 - vz2) / h;
        tzz_z[ off1 + i*Npml ] *= c1;
		for (size_t j=1; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]   + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + i*Npml+j-1] + vz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[j] );
			c2 = ( 2. - dt*d[j] ) / (2.*dt);
			txx_z[ off1 + i*Npml+j ] *= c2;
			txx_z[ off1 + i*Npml+j ] += la[(nx-1)*nz] * (vz1 - vz2) / h;
			txx_z[ off1 + i*Npml+j ] *= c1;
			tzz_z[ off1 + i*Npml+j ] *= c2;
			tzz_z[ off1 + i*Npml+j ] += l2m[(nx-1)*nz] * (vz1 - vz2) / h;
			tzz_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	
	// txx & tzz - right strip
	//
	off1 += Npml*Npml;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<nz; ++j) {
			vx1 = vx_x[off1 + (i+1)*nz+j] + vx_z[off1 + (i+1)*nz+j];
			vx2 = vx_x[off1 +     i*nz+j] + vx_z[off1 +     i*nz+j];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-i] );
			c2 = ( 2. - dt*dH[Npml-1-i] ) / (2.*dt);
			txx_x[ off1 + i*nz+j ] *= c2;
			txx_x[ off1 + i*nz+j ] += l2m[(nx-1)*nz+j] * (vx1 - vx2) / h;
			txx_x[ off1 + i*nz+j ] *= c1;
			tzz_x[ off1 + i*nz+j ] *= c2;
			tzz_x[ off1 + i*nz+j ] += la[(nx-1)*nz+j] * (vx1 - vx2) / h;
			tzz_x[ off1 + i*nz+j ] *= c1;
		}
	}
    for (size_t j=0; j<nz; ++j) {
        vx1 = 0.;
        vx2 = vx_x[off1 + (Npml-1)*nz+j] + vx_z[off1 + (Npml-1)*nz+j];
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        txx_x[ off1 + (Npml-1)*nz+j ] *= c2;
        txx_x[ off1 + (Npml-1)*nz+j ] += l2m[(nx-1)*nz+j] * (vx1 - vx2) / h;
        txx_x[ off1 + (Npml-1)*nz+j ] *= c1;
        tzz_x[ off1 + (Npml-1)*nz+j ] *= c2;
        tzz_x[ off1 + (Npml-1)*nz+j ] += la[(nx-1)*nz+j] * (vx1 - vx2) / h;
        tzz_x[ off1 + (Npml-1)*nz+j ] *= c1;
    }
	off2 = 2*Npml*Npml + Npml*nz + 2*Npml*nx;
	for (size_t i=0; i<Npml; ++i) {
		vz1 = vz_x[off1 + i*nz]          + vz_z[off1 + i*nz];
		vz2 = vz_x[off2 + i*Npml+Npml-1] + vz_z[off2 + i*Npml+Npml-1];
		txx_z[ off1 + i*nz ] += dt * la[(nx-1)*nz] * (vz1 - vz2) / h;
		tzz_z[ off1 + i*nz ] += dt * l2m[(nx-1)*nz] * (vz1 - vz2) / h;
		for (size_t j=1; j<nz; ++j) {
			vz1 = vz_x[off1 + i*nz+j]   + vz_z[off1 + i*nz+j];
			vz2 = vz_x[off1 + i*nz+j-1] + vz_z[off1 + i*nz+j-1];
			txx_z[ off1 + i*nz+j ] += dt * la[(nx-1)*nz+j] * (vz1 - vz2) / h;
			tzz_z[ off1 + i*nz+j ] += dt * l2m[(nx-1)*nz+j] * (vz1 - vz2) / h;
		}
	}
	
	// txx & tzz - lower right corner
	//
	off1 += Npml*nz;
	for (size_t i=0; i<Npml-1; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vx1 = vx_x[off1 + (i+1)*Npml+j] + vx_z[off1 + (i+1)*Npml+j];
			vx2 = vx_x[off1 +     i*Npml+j] + vx_z[off1 +     i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-i] );
			c2 = ( 2. - dt*dH[Npml-1-i] ) / (2.*dt);
			txx_x[ off1 + i*Npml+j ] *= c2;
			txx_x[ off1 + i*Npml+j ] += l2m[nx*nx-1] * (vx1 - vx2) / h;
			txx_x[ off1 + i*Npml+j ] *= c1;
			tzz_x[ off1 + i*Npml+j ] *= c2;
			tzz_x[ off1 + i*Npml+j ] += la[nx*nx-1] * (vx1 - vx2) / h;
			tzz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
    for (size_t j=0; j<Npml; ++j) {
        vx1 = 0.;
        vx2 = vx_x[off1 + (Npml-1)*Npml+j] + vx_z[off1 + (Npml-1)*Npml+j];
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        txx_x[ off1 + (Npml-1)*Npml+j ] *= c2;
        txx_x[ off1 + (Npml-1)*Npml+j ] += l2m[nx*nx-1] * (vx1 - vx2) / h;
        txx_x[ off1 + (Npml-1)*Npml+j ] *= c1;
        tzz_x[ off1 + (Npml-1)*Npml+j ] *= c2;
        tzz_x[ off1 + (Npml-1)*Npml+j ] += la[nx*nx-1] * (vx1 - vx2) / h;
        tzz_x[ off1 + (Npml-1)*Npml+j ] *= c1;
    }
	off2 = 3*Npml*Npml + Npml*nz + 2*Npml*nx;
	for (size_t i=0; i<Npml; ++i) {
		vz1 = vz_x[off1 + i*Npml]    + vz_z[off1 + i*Npml];
		vz2 = vz_x[off2 + i*nz+nz-1] + vz_z[off2 + i*nz+nz-1];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		txx_z[ off1 + i*Npml ] *= c2;
		txx_z[ off1 + i*Npml ] += la[nx*nx-1] * (vz1 - vz2) / h;
		txx_z[ off1 + i*Npml ] *= c1;
		tzz_z[ off1 + i*Npml ] *= c2;
		tzz_z[ off1 + i*Npml ] += l2m[nx*nx-1] * (vz1 - vz2) / h;
		tzz_z[ off1 + i*Npml ] *= c1;
		for (size_t j=1; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]   + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + i*Npml+j-1] + vz_z[off1 + i*Npml+j-1];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-j] );
			c2 = ( 2. - dt*d[Npml-1-j] ) / (2.*dt);
			txx_z[ off1 + i*Npml+j ] *= c2;
			txx_z[ off1 + i*Npml+j ] += la[nx*nx-1] * (vz1 - vz2) / h;
			txx_z[ off1 + i*Npml+j ] *= c1;
			tzz_z[ off1 + i*Npml+j ] *= c2;
			tzz_z[ off1 + i*Npml+j ] += l2m[nx*nx-1] * (vz1 - vz2) / h;
			tzz_z[ off1 + i*Npml+j ] *= c1;
		}
	}
	
    
    
    
    // txz - upper left corner
	//
	off1 = 0;
    for (size_t j=0; j<Npml; ++j) {
        vz1 = vz_x[off1 + j] + vz_z[off1 + j];
        vz2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        txz_x[ off1 + j ] *= c2;
        txz_x[ off1 + j ] += mu[0] * (vz1 - vz2) / h;
        txz_x[ off1 + j ] *= c1;
    }
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]     + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + (i-1)*Npml+j] + vz_z[off1 + (i-1)*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*d[i] );
			c2 = ( 2. - dt*d[i] ) / (2.*dt);
			txz_x[ off1 + i*Npml+j ] *= c2;
			txz_x[ off1 + i*Npml+j ] += mu[0] * (vz1 - vz2) / h;
			txz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
	off2 = Npml*Npml;
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			vx1 = vx_x[off1 + i*Npml+j+1] + vx_z[off1 + i*Npml+j+1];
			vx2 = vx_x[off1 + i*Npml+j]   + vx_z[off1 + i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dh[j] );
			c2 = ( 2. - dt*dh[j] ) / (2.*dt);
			txz_z[ off1 + i*Npml+j ] *= c2;
			txz_z[ off1 + i*Npml+j ] += mu[0] * (vx1 - vx2) / h;
			txz_z[ off1 + i*Npml+j ] *= c1;
		}
		vx1 = vx_x[off2 + i*nz]          + vx_z[off2 + i*nz];
		vx2 = vx_x[off1 + i*Npml+Npml-1] + vx_z[off1 + i*Npml+Npml-1];
		c1 = (2.*dt) / ( 2. + dt*dh[Npml-1] );
		c2 = ( 2. - dt*dh[Npml-1] ) / (2.*dt);
		txz_z[ off1 + i*Npml+Npml-1 ] *= c2;
		txz_z[ off1 + i*Npml+Npml-1 ] += mu[0] * (vx1 - vx2) / h;
		txz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}
	
	// txz - left strip
	//
	off1 += Npml*Npml;
    for (size_t j=0; j<nz; ++j) {
        vz1 = vz_x[off1 + j] + vz_z[off1 + j];
        vz2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        txz_x[ off1 + j ] *= c2;
        txz_x[ off1 + j ] += mu[j] * (vz1 - vz2) / h;
        txz_x[ off1 + j ] *= c1;
    }
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<nz; ++j) {
			vz1 = vz_x[off1 + i*nz+j]     + vz_z[off1 + i*nz+j];
			vz2 = vz_x[off1 + (i-1)*nz+j] + vz_z[off1 + (i-1)*nz+j];
			c1 = (2.*dt) / ( 2. + dt*d[i] );
			c2 = ( 2. - dt*d[i] ) / (2.*dt);
			txz_x[ off1 + i*nz+j ] *= c2;
			txz_x[ off1 + i*nz+j ] += mu[j] * (vz1 - vz2) / h;
			txz_x[ off1 + i*nz+j ] *= c1;
		}
	}
	off2 = Npml*Npml + Npml*nz;
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<nz-1; ++j) {
			vx1 = vx_x[off1 + i*nz+j+1] + vx_z[off1 + i*nz+j+1];
			vx2 = vx_x[off1 + i*nz+j]   + vx_z[off1 + i*nz+j];
			txz_z[ off1 + i*nz+j ] += dt * mu[j] * (vx1 - vx2) / h;
		}
		vx1 = vx_x[off2 + i*Npml]    + vx_z[off2 + i*Npml];
		vx2 = vx_x[off1 + i*nz+nz-1] + vx_z[off1 + i*nz+nz-1];
		txz_z[ off1 + i*nz+nz-1 ] += dt * mu[nz-1] * (vx1 - vx2) / h;
	}
	
	// txz - lower left corner
	//
	off1 += Npml*nz;
    for (size_t j=0; j<Npml; ++j) {
        vz1 = vz_x[off1 + j] + vz_z[off1 + j];
        vz2 = 0.;
        c1 = (2.*dt) / ( 2. + dt*d[0] );
        c2 = ( 2. - dt*d[0] ) / (2.*dt);
        txz_x[ off1 + j ] *= c2;
        txz_x[ off1 + j ] += mu[nz-1] * (vz1 - vz2) / h;
        txz_x[ off1 + j ] *= c1;
    }
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]     + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + (i-1)*Npml+j] + vz_z[off1 + (i-1)*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*d[i] );
			c2 = ( 2. - dt*d[i] ) / (2.*dt);
			txz_x[ off1 + i*Npml+j ] *= c2;
			txz_x[ off1 + i*Npml+j ] += mu[nz-1] * (vz1 - vz2) / h;
			txz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			vx1 = vx_x[off1 + i*Npml+j+1] + vx_z[off1 + i*Npml+j+1];
			vx2 = vx_x[off1 + i*Npml+j]   + vx_z[off1 + i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-j] );
			c2 = ( 2. - dt*dH[Npml-1-j] ) / (2.*dt);
			txz_z[ off1 + i*Npml+j ] *= c2;
			txz_z[ off1 + i*Npml+j ] += mu[nz-1] * (vx1 - vx2) / h;
			txz_z[ off1 + i*Npml+j ] *= c1;
		}
        vx1 = 0.;
        vx2 = vx_x[off1 + i*Npml+Npml-1]   + vx_z[off1 + i*Npml+Npml-1];
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        txz_z[ off1 + i*Npml+Npml-1 ] *= c2;
        txz_z[ off1 + i*Npml+Npml-1 ] += mu[nz-1] * (vx1 - vx2) / h;
        txz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}
	
	// txz - upper strip
	//
	off1 += Npml*Npml;
	off2 = 0;
	for (size_t j=0; j<Npml; ++j) {
		vz1 = vz_x[off1 + j]               + vz_z[off1 + j];
		vz2 = vz_x[off2 + (Npml-1)*Npml+j] + vz_z[off2 + (Npml-1)*Npml+j];
		txz_x[ off1 + j ] += dt * mu[0] * (vz1 - vz2) / h;
	}
	for (size_t i=1; i<nx; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]     + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + (i-1)*Npml+j] + vz_z[off1 + (i-1)*Npml+j];
			txz_x[ off1 + i*Npml+j ] += dt * mu[i*nz] * (vz1 - vz2) / h;
		}
	}
	for (size_t i=0; i<nx; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			vx1 = vx_x[off1 + i*Npml+j+1] + vx_z[off1 + i*Npml+j+1];
			vx2 = vx_x[off1 + i*Npml+j]   + vx_z[off1 + i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dh[j] );
			c2 = ( 2. - dt*dh[j] ) / (2.*dt);
			txz_z[ off1 + i*Npml+j ] *= c2;
			txz_z[ off1 + i*Npml+j ] += mu[i*nz] * (vx1 - vx2) / h;
			txz_z[ off1 + i*Npml+j ] *= c1;
		}
		vx1 = vx[i*nz];
		vx2 = vx_x[off1 + i*Npml+Npml-1]   + vx_z[off1 + i*Npml+Npml-1];
		c1 = (2.*dt) / ( 2. + dt*dh[Npml-1] );
		c2 = ( 2. - dt*dh[Npml-1] ) / (2.*dt);
		txz_z[ off1 + i*Npml+Npml-1 ] *= c2;
		txz_z[ off1 + i*Npml+Npml-1 ] += mu[i*nz] * (vx1 - vx2) / h;
		txz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}
	
	// txz - lower strip
	//
	off1 += Npml*nx;
	off2 = Npml*Npml + Npml*nz;
	for (size_t j=0; j<Npml; ++j) {
		vz1 = vz_x[off1 + j]               + vz_z[off1 + j];
		vz2 = vz_x[off2 + (Npml-1)*Npml+j] + vz_z[off2 + (Npml-1)*Npml+j];
		txz_x[ off1 + j ] += dt * mu[nz-1] * (vz1 - vz2) / h;
	}
	for (size_t i=1; i<nx; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]     + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + (i-1)*Npml+j] + vz_z[off1 + (i-1)*Npml+j];
			txz_x[ off1 + i*Npml+j ] += dt * mu[i*nz+nz-1] * (vz1 - vz2) / h;
		}
	}
	for (size_t i=0; i<nx; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			vx1 = vx_x[off1 + i*Npml+j+1] + vx_z[off1 + i*Npml+j+1];
			vx2 = vx_x[off1 + i*Npml+j]   + vx_z[off1 + i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-j] );
			c2 = ( 2. - dt*dH[Npml-1-j] ) / (2.*dt);
			txz_z[ off1 + i*Npml+j ] *= c2;
			txz_z[ off1 + i*Npml+j ] += mu[i*nz+nz-1] * (vx1 - vx2) / h;
			txz_z[ off1 + i*Npml+j ] *= c1;
		}
        vx1 = 0.;
        vx2 = vx_x[off1 + i*Npml+Npml-1] + vx_z[off1 + i*Npml+Npml-1];
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        txz_z[ off1 + i*Npml+Npml-1 ] *= c2;
        txz_z[ off1 + i*Npml+Npml-1 ] += mu[i*nz+nz-1] * (vx1 - vx2) / h;
        txz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}
	
	// txz - upper right corner
	//
	off1 += Npml*nx;
	off2 = 2*Npml*Npml + Npml*nz;
	for (size_t j=0; j<Npml; ++j) {
		vz1 = vz_x[off1 + j]             + vz_z[off1 + j];
		vz2 = vz_x[off2 + (nx-1)*Npml+j] + vz_z[off2 + (nx-1)*Npml+j];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		txz_x[ off1 + j ] *= c2;
		txz_x[ off1 + j ] += mu[(nx-1)*nz] * (vz1 - vz2) / h;
		txz_x[ off1 + j ] *= c1;
	}
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]     + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + (i-1)*Npml+j] + vz_z[off1 + (i-1)*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-i] );
			c2 = ( 2. - dt*d[Npml-1-i] ) / (2.*dt);
			txz_x[ off1 + i*Npml+j ] *= c2;
			txz_x[ off1 + i*Npml+j ] += mu[(nx-1)*nz] * (vz1 - vz2) / h;
			txz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
	off2 = 3*Npml*Npml + Npml*nz + 2*Npml*nx;
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			vx1 = vx_x[off1 + i*Npml+j+1] + vx_z[off1 + i*Npml+j+1];
			vx2 = vx_x[off1 + i*Npml+j]   + vx_z[off1 + i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dh[j] );
			c2 = ( 2. - dt*dh[j] ) / (2.*dt);
			txz_z[ off1 + i*Npml+j ] *= c2;
			txz_z[ off1 + i*Npml+j ] += mu[(nx-1)*nz] * (vx1 - vx2) / h;
			txz_z[ off1 + i*Npml+j ] *= c1;
		}
		vx1 = vx_x[off2 + i*nz]          + vx_z[off2 + i*nz];
		vx2 = vx_x[off1 + i*Npml+Npml-1] + vx_z[off1 + i*Npml+Npml-1];
		c1 = (2.*dt) / ( 2. + dt*dh[Npml-1] );
		c2 = ( 2. - dt*dh[Npml-1] ) / (2.*dt);
		txz_z[ off1 + i*Npml+Npml-1 ] *= c2;
		txz_z[ off1 + i*Npml+Npml-1 ] += mu[(nx-1)*nz] * (vx1 - vx2) / h;
		txz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}
	
	// txz - right strip
	//
	off1 += Npml*Npml;
	for (size_t j=0; j<nz; ++j) {
		vz1 = vz_x[off1 + j]   + vz_z[off1 + j];
		vz2 = vz[(nx-1)*nz+j];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		txz_x[ off1 + j ] *= c2;
		txz_x[ off1 + j ] += mu[(nx-1)*nz+j] * (vz1 - vz2) / h;
		txz_x[ off1 + j ] *= c1;
	}
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<nz; ++j) {
			vz1 = vz_x[off1 + i*nz+j]     + vz_z[off1 + i*nz+j];
			vz2 = vz_x[off1 + (i-1)*nz+j] + vz_z[off1 + (i-1)*nz+j];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-i] );
			c2 = ( 2. - dt*d[Npml-1-i] ) / (2.*dt);
			txz_x[ off1 + i*nz+j ] *= c2;
			txz_x[ off1 + i*nz+j ] += mu[(nx-1)*nz+j] * (vz1 - vz2) / h;
			txz_x[ off1 + i*nz+j ] *= c1;
		}
	}
	off2 = 3*Npml*Npml + 2*Npml*nz + 2*Npml*nx;
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<nz-1; ++j) {
			vx1 = vx_x[off1 + i*nz+j+1] + vx_z[off1 + i*nz+j+1];
			vx2 = vx_x[off1 + i*nz+j]   + vx_z[off1 + i*nz+j];
			txz_z[ off1 + i*nz+j ] += dt * mu[(nx-1)*nz+j] * (vx1 - vx2) / h;
		}
		vx1 = vx_x[off2 + i*Npml]    + vx_z[off2 + i*Npml];
		vx2 = vx_x[off1 + i*nz+nz-1] + vx_z[off1 + i*nz+nz-1];
		txz_z[ off1 + i*nz+nz-1 ] += dt * mu[nx*nz-1] * (vx1 - vx2) / h;
	}
	
	// txz - lower right corner
	//
	off1 += Npml*nz;
	off2 = 2*Npml*Npml + Npml*nz + Npml*nx;
	for (size_t j=0; j<Npml; ++j) {
		vz1 = vz_x[off1 + j]             + vz_z[off1 + j];
		vz2 = vz_x[off2 + (nx-1)*Npml+j] + vz_z[off2 + (nx-1)*Npml+j];
		c1 = (2.*dt) / ( 2. + dt*d[Npml-1] );
		c2 = ( 2. - dt*d[Npml-1] ) / (2.*dt);
		txz_x[ off1 + j ] *= c2;
		txz_x[ off1 + j ] += mu[nx*nz-1] * (vz1 - vz2) / h;
		txz_x[ off1 + j ] *= c1;
	}
	for (size_t i=1; i<Npml; ++i) {
		for (size_t j=0; j<Npml; ++j) {
			vz1 = vz_x[off1 + i*Npml+j]     + vz_z[off1 + i*Npml+j];
			vz2 = vz_x[off1 + (i-1)*Npml+j] + vz_z[off1 + (i-1)*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*d[Npml-1-i] );
			c2 = ( 2. - dt*d[Npml-1-i] ) / (2.*dt);
			txz_x[ off1 + i*Npml+j ] *= c2;
			txz_x[ off1 + i*Npml+j ] += mu[nx*nz-1] * (vz1 - vz2) / h;
			txz_x[ off1 + i*Npml+j ] *= c1;
		}
	}
	for (size_t i=0; i<Npml; ++i) {
		for (size_t j=0; j<Npml-1; ++j) {
			vx1 = vx_x[off1 + i*Npml+j+1] + vx_z[off1 + i*Npml+j+1];
			vx2 = vx_x[off1 + i*Npml+j]   + vx_z[off1 + i*Npml+j];
			c1 = (2.*dt) / ( 2. + dt*dH[Npml-1-j] );
			c2 = ( 2. - dt*dH[Npml-1-j] ) / (2.*dt);
			txz_z[ off1 + i*Npml+j ] *= c2;
			txz_z[ off1 + i*Npml+j ] += mu[nx*nz-1] * (vx1 - vx2) / h;
			txz_z[ off1 + i*Npml+j ] *= c1;
		}
        vx1 = 0.;
        vx2 = vx_x[off1 + i*Npml+Npml-1] + vx_z[off1 + i*Npml+Npml-1];
        c1 = (2.*dt) / ( 2. + dt*dH[0] );
        c2 = ( 2. - dt*dH[0] ) / (2.*dt);
        txz_z[ off1 + i*Npml+Npml-1 ] *= c2;
        txz_z[ off1 + i*Npml+Npml-1 ] += mu[nx*nz-1] * (vx1 - vx2) / h;
        txz_z[ off1 + i*Npml+Npml-1 ] *= c1;
	}
}	

void alloc_cpml(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
				struct mem_pml mem[3], const struct grid *g) {
    if ( NULL == ( fp1->k_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->k_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kh_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kh_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->b_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->c_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->bh_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->ch_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->b_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->c_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->bh_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->ch_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->b_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->c_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->bh_x = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->ch_x = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->b_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->c_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->bh_z = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->ch_z = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->b_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->c_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->bh_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->ch_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->b_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->c_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->bh_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->ch_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->bH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->cH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->bH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->cH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->bH_x = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->cH_x = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->bH_z = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->cH_z = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->bH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->cH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->bH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->cH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    for ( size_t n=0; n<3; ++n ) {
		if ( NULL == ( mem[n].dx_txx = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dx_p   = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dx_txz = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dx_vx  = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dx_qx  = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dx_vz  = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_txz = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_tzz = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_p   = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_vz  = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_qz  = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_vx  = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	}
	fp4->k_x  = fp23->k_x  = fp1->k_x;
	fp4->k_z  = fp23->k_z  = fp1->k_z;
	fp4->kh_x = fp23->kh_x = fp1->kh_x;
	fp4->kh_z = fp23->kh_z = fp1->kh_z;
	fp4->kH_x = fp23->kH_x = fp1->kH_x;
	fp4->kH_z = fp23->kH_z = fp1->kH_z;
	
	for ( size_t n=0; n<(2*g->ab.np+1)*g->nz2; ++n ) {
        mem[0].dx_txx[n] = mem[0].dx_p[n]  = mem[0].dx_txz[n] = 0.0;
        mem[0].dx_vx[n]  = mem[0].dx_qx[n] = mem[0].dx_vz[n] = 0.0;
        mem[1].dx_txx[n] = mem[1].dx_p[n]  = mem[1].dx_txz[n] = 0.0;
        mem[1].dx_vx[n]  = mem[1].dx_qx[n] = mem[1].dx_vz[n] = 0.0;
        mem[2].dx_txx[n] = mem[2].dx_p[n]  = mem[2].dx_txz[n] = 0.0;
        mem[2].dx_vx[n]  = mem[2].dx_qx[n] = mem[2].dx_vz[n] = 0.0;
    }        
    for ( size_t n=0; n<(2*g->ab.np+1)*g->nx2; ++n ) {
        mem[0].dz_txz[n] = mem[0].dz_tzz[n] = mem[0].dz_p[n] = 0.0;
        mem[0].dz_vz[n]  = mem[0].dz_qz[n]  = mem[0].dz_vx[n] = 0.0; 
        mem[1].dz_txz[n] = mem[1].dz_tzz[n] = mem[1].dz_p[n] = 0.0;
        mem[1].dz_vz[n]  = mem[1].dz_qz[n]  = mem[1].dz_vx[n] = 0.0; 
        mem[2].dz_txz[n] = mem[2].dz_tzz[n] = mem[2].dz_p[n] = 0.0;
        mem[2].dz_vz[n]  = mem[2].dz_qz[n]  = mem[2].dz_vx[n] = 0.0; 
    }
	
}

void free_cpml(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
				struct mem_pml mem[3]) {
    free ( fp1->k_x );
    free ( fp1->k_z );
    free ( fp1->kh_x );
    free ( fp1->kh_z );
    free ( fp1->b_x );
    free ( fp1->c_x );
    free ( fp1->bh_x );
    free ( fp1->ch_x );
    free ( fp1->b_z );
    free ( fp1->c_z );
    free ( fp1->bh_z );
    free ( fp1->ch_z );
    free ( fp23->b_x );
    free ( fp23->c_x );
    free ( fp23->bh_x );
    free ( fp23->ch_x );
    free ( fp23->b_z );
    free ( fp23->c_z );
    free ( fp23->bh_z );
    free ( fp23->ch_z );
    free ( fp4->b_x );
    free ( fp4->c_x );
    free ( fp4->bh_x );
    free ( fp4->ch_x );
    free ( fp4->b_z );
    free ( fp4->c_z );
    free ( fp4->bh_z );
    free ( fp4->ch_z );
	free ( fp1->kH_x );
    free ( fp1->kH_z );
	free ( fp1->bH_x );
	free ( fp1->cH_x );
	free ( fp1->bH_z );
	free ( fp1->cH_z );
	free ( fp23->bH_x );
	free ( fp23->cH_x );
	free ( fp23->bH_z );
	free ( fp23->cH_z );
	free ( fp4->bH_x );
	free ( fp4->cH_x );
	free ( fp4->bH_z );
	free ( fp4->cH_z );
	for ( size_t n=0; n<3; ++n ) {
		free ( mem[n].dx_txx );
		free ( mem[n].dx_p );
		free ( mem[n].dx_txz );
		free ( mem[n].dx_vx );
		free ( mem[n].dx_qx );
		free ( mem[n].dx_vz );
		free ( mem[n].dz_txz );
		free ( mem[n].dz_tzz );
		free ( mem[n].dz_p );
		free ( mem[n].dz_vz );
		free ( mem[n].dz_qz );
		free ( mem[n].dz_vx );
	}
}

void alloc_cpml_ve(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
                   struct mem_pml mem[3], const struct grid *g) {
    if ( NULL == ( fp1->k_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->k_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kh_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kh_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->b_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->c_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->bh_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->ch_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->b_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->c_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->bh_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->ch_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->b_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->c_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->bh_x = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->ch_x = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->b_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->c_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->bh_z = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->ch_z = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->b_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->c_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->bh_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->ch_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->b_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->c_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->bh_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->ch_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->bH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->cH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->bH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->cH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->bH_x = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->cH_x = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->bH_z = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->cH_z = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->bH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->cH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->bH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->cH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    for ( size_t n=0; n<3; ++n ) {
		if ( NULL == ( mem[n].dx_txx = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dx_txz = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dx_vx  = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dx_vz  = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_txz = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_tzz = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_vz  = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_vx  = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	}
	fp4->k_x  = fp23->k_x  = fp1->k_x;
	fp4->k_z  = fp23->k_z  = fp1->k_z;
	fp4->kh_x = fp23->kh_x = fp1->kh_x;
	fp4->kh_z = fp23->kh_z = fp1->kh_z;
	fp4->kH_x = fp23->kH_x = fp1->kH_x;
	fp4->kH_z = fp23->kH_z = fp1->kH_z;
	
	for ( size_t n=0; n<(2*g->ab.np+1)*g->nz2; ++n ) {
        mem[0].dx_txx[n] = mem[0].dx_txz[n] = 0.0;
        mem[0].dx_vx[n]  = mem[0].dx_vz[n]  = 0.0;
        mem[1].dx_txx[n] = mem[1].dx_txz[n] = 0.0;
        mem[1].dx_vx[n]  = mem[1].dx_vz[n]  = 0.0;
        mem[2].dx_txx[n] = mem[2].dx_txz[n] = 0.0;
        mem[2].dx_vx[n]  = mem[2].dx_vz[n]  = 0.0;
    }
    for ( size_t n=0; n<(2*g->ab.np+1)*g->nx2; ++n ) {
        mem[0].dz_txz[n] = mem[0].dz_tzz[n] = 0.0;
        mem[0].dz_vz[n]  = mem[0].dz_vx[n]  = 0.0;
        mem[1].dz_txz[n] = mem[1].dz_tzz[n] = 0.0;
        mem[1].dz_vz[n]  = mem[1].dz_vx[n]  = 0.0;
        mem[2].dz_txz[n] = mem[2].dz_tzz[n] = 0.0;
        mem[2].dz_vz[n]  = mem[2].dz_vx[n]  = 0.0;
    }
	
}

void free_cpml_ve(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
                  struct mem_pml mem[3]) {
    free ( fp1->k_x );
    free ( fp1->k_z );
    free ( fp1->kh_x );
    free ( fp1->kh_z );
    free ( fp1->b_x );
    free ( fp1->c_x );
    free ( fp1->bh_x );
    free ( fp1->ch_x );
    free ( fp1->b_z );
    free ( fp1->c_z );
    free ( fp1->bh_z );
    free ( fp1->ch_z );
    free ( fp23->b_x );
    free ( fp23->c_x );
    free ( fp23->bh_x );
    free ( fp23->ch_x );
    free ( fp23->b_z );
    free ( fp23->c_z );
    free ( fp23->bh_z );
    free ( fp23->ch_z );
    free ( fp4->b_x );
    free ( fp4->c_x );
    free ( fp4->bh_x );
    free ( fp4->ch_x );
    free ( fp4->b_z );
    free ( fp4->c_z );
    free ( fp4->bh_z );
    free ( fp4->ch_z );
	free ( fp1->kH_x );
    free ( fp1->kH_z );
	free ( fp1->bH_x );
	free ( fp1->cH_x );
	free ( fp1->bH_z );
	free ( fp1->cH_z );
	free ( fp23->bH_x );
	free ( fp23->cH_x );
	free ( fp23->bH_z );
	free ( fp23->cH_z );
	free ( fp4->bH_x );
	free ( fp4->cH_x );
	free ( fp4->bH_z );
	free ( fp4->cH_z );
	for ( size_t n=0; n<3; ++n ) {
		free ( mem[n].dx_txx );
		free ( mem[n].dx_txz );
		free ( mem[n].dx_vx );
		free ( mem[n].dx_vz );
		free ( mem[n].dz_txz );
		free ( mem[n].dz_tzz );
		free ( mem[n].dz_vz );
		free ( mem[n].dz_vx );
	}
}

void alloc_cpml_ve_sh(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
					  struct mem_pml_sh mem[3], const struct grid *g) {
    if ( NULL == ( fp1->k_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->k_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kh_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kh_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->b_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->c_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->bh_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->ch_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->b_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->c_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->bh_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->ch_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->b_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->c_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->bh_x = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->ch_x = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->b_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->c_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->bh_z = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp23->ch_z = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->b_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->c_x   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->bh_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->ch_x  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->b_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->c_z   = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->bh_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp4->ch_z  = (double *) malloc( g->ab.np * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    if ( NULL == ( fp1->kH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->bH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->cH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->bH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp1->cH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->bH_x = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->cH_x = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->bH_z = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp23->cH_z = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->bH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->cH_x  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->bH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	if ( NULL == ( fp4->cH_z  = (double *) malloc( (g->ab.np+1) * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    for ( size_t n=0; n<3; ++n ) {
		if ( NULL == ( mem[n].dx_txy = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dx_vy  = (double *) malloc( (2*g->ab.np+1)*g->nz2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_tyz = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( mem[n].dz_vy  = (double *) malloc( (2*g->ab.np+1)*g->nx2 * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	}
	fp4->k_x  = fp23->k_x  = fp1->k_x;
	fp4->k_z  = fp23->k_z  = fp1->k_z;
	fp4->kh_x = fp23->kh_x = fp1->kh_x;
	fp4->kh_z = fp23->kh_z = fp1->kh_z;
	fp4->kH_x = fp23->kH_x = fp1->kH_x;
	fp4->kH_z = fp23->kH_z = fp1->kH_z;
	
	for ( size_t n=0; n<(2*g->ab.np+1)*g->nz2; ++n ) {
        mem[0].dx_txy[n] = mem[0].dx_vy[n]  = 0.0;
        mem[1].dx_txy[n] = mem[1].dx_vy[n]  = 0.0;
        mem[2].dx_txy[n] = mem[2].dx_vy[n]  = 0.0;
    }
    for ( size_t n=0; n<(2*g->ab.np+1)*g->nx2; ++n ) {
        mem[0].dz_tyz[n] = mem[0].dz_vy[n]  = 0.0;
        mem[1].dz_tyz[n] = mem[1].dz_vy[n]  = 0.0;
        mem[2].dz_tyz[n] = mem[2].dz_vy[n]  = 0.0;
    }
	
}

void free_cpml_ve_sh(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
					 struct mem_pml_sh mem[3]) {
    free ( fp1->k_x );
    free ( fp1->k_z );
    free ( fp1->kh_x );
    free ( fp1->kh_z );
    free ( fp1->b_x );
    free ( fp1->c_x );
    free ( fp1->bh_x );
    free ( fp1->ch_x );
    free ( fp1->b_z );
    free ( fp1->c_z );
    free ( fp1->bh_z );
    free ( fp1->ch_z );
    free ( fp23->b_x );
    free ( fp23->c_x );
    free ( fp23->bh_x );
    free ( fp23->ch_x );
    free ( fp23->b_z );
    free ( fp23->c_z );
    free ( fp23->bh_z );
    free ( fp23->ch_z );
    free ( fp4->b_x );
    free ( fp4->c_x );
    free ( fp4->bh_x );
    free ( fp4->ch_x );
    free ( fp4->b_z );
    free ( fp4->c_z );
    free ( fp4->bh_z );
    free ( fp4->ch_z );
	free ( fp1->kH_x );
    free ( fp1->kH_z );
	free ( fp1->bH_x );
	free ( fp1->cH_x );
	free ( fp1->bH_z );
	free ( fp1->cH_z );
	free ( fp23->bH_x );
	free ( fp23->cH_x );
	free ( fp23->bH_z );
	free ( fp23->cH_z );
	free ( fp4->bH_x );
	free ( fp4->cH_x );
	free ( fp4->bH_z );
	free ( fp4->cH_z );
	for ( size_t n=0; n<3; ++n ) {
		free ( mem[n].dx_txy );
		free ( mem[n].dx_vy );
		free ( mem[n].dz_tyz );
		free ( mem[n].dz_vy );
	}
}

void compute_cpml(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
				  const struct grid *g, const double dt, const double alpha_max,
				  const short iwipe) {

//	const double pi = 4.0*atan(1.0);
    double Lx = g->ab.np*g->dx;
    double Lz = g->ab.np*g->dz;
    double d0x = -(g->ab.pmlOrder + 1.) * g->ab.vmax * log(g->ab.Rc)/(2.*Lx);
    double d0z = -(g->ab.pmlOrder + 1.) * g->ab.vmax * log(g->ab.Rc)/(2.*Lz);
    double da = alpha_max/g->ab.np;
	double theta = 0.5;  // semi-implicit scheme
    
	if ( iwipe ) {
		for ( size_t n=0; n<g->ab.np; ++n ) {
			
			fp1->k_x[n] = 1.+(g->ab.kappa_max-1.) * pow(((g->ab.np-n)*g->dx/Lx), g->ab.pmlOrder);
			fp1->k_z[n] = fp1->k_x[n];
			fp1->kh_x[n] = 1.+(g->ab.kappa_max-1.) * pow(((g->ab.np-n-0.5)*g->dx/Lx), g->ab.pmlOrder);
			fp1->kh_z[n] = fp1->kh_x[n];
			
			double d_x = d0x * pow(((g->ab.np-n)*g->dx/Lx), g->ab.pmlOrder);
			double alpha_x = n*da;
			double dh_x = d0x * pow(((g->ab.np-n-0.5)*g->dx/Lx), g->ab.pmlOrder);
			double alphah_x = (n+0.5)*da;
			
			fp1->b_x[n] = 1.0;
			fp1->c_x[n] = 0.0;
			fp1->bh_x[n] = 1.0;
			fp1->ch_x[n] = 0.0;
			
			fp23->b_x[n] = ( fp1->k_x[n] - (1.-theta)*dt*0.5*(d_x+fp1->k_x[n]*alpha_x) ) /
			( fp1->k_x[n] +      theta*dt*0.5*(d_x+fp1->k_x[n]*alpha_x) );
			fp23->bh_x[n] = ( fp1->kh_x[n] - (1.-theta)*dt*0.5*(dh_x+fp1->kh_x[n]*alphah_x) ) /
			( fp1->kh_x[n] +      theta*dt*0.5*(dh_x+fp1->kh_x[n]*alphah_x) );
			fp23->c_x[n] = -d_x * dt*0.5 /
			( fp1->k_x[n]*(fp1->k_x[n]+theta*dt*0.5*(d_x+fp1->k_x[n]*alpha_x)));
			fp23->ch_x[n] = -dh_x * dt*0.5 /
			( fp1->kh_x[n]*(fp1->kh_x[n]+theta*dt*0.5*(dh_x+fp1->kh_x[n]*alphah_x)));
			
			fp4->b_x[n] = ( fp1->k_x[n] - (1.-theta)*dt*(d_x+fp1->k_x[n]*alpha_x) ) /
			( fp1->k_x[n] +      theta*dt*(d_x+fp1->k_x[n]*alpha_x) );
			fp4->bh_x[n] = ( fp1->kh_x[n] - (1.-theta)*dt*(dh_x+fp1->kh_x[n]*alphah_x) ) /
			( fp1->kh_x[n] +      theta*dt*(dh_x+fp1->kh_x[n]*alphah_x) );
			fp4->c_x[n] = -d_x * dt /
			( fp1->k_x[n]*(fp1->k_x[n]+theta*dt*(d_x+fp1->k_x[n]*alpha_x)));
			fp4->ch_x[n] = -dh_x * dt /
			( fp1->kh_x[n]*(fp1->kh_x[n]+theta*dt*(dh_x+fp1->kh_x[n]*alphah_x)));
			
			double d_z = d0z * pow(((g->ab.np-n)*g->dz/Lz), g->ab.pmlOrder);
			double alpha_z = n*da;
			double dh_z = d0z * pow(((g->ab.np-n-0.5)*g->dz/Lz), g->ab.pmlOrder);
			double alphah_z = (n+0.5)*da;
			
			fp1->b_z[n] = 1.0;
			fp1->c_z[n] = 0.0;
			fp1->bh_z[n] = 1.0;
			fp1->ch_z[n] = 0.0;
			
			fp23->b_z[n] = ( fp1->k_z[n] - (1.-theta)*dt*0.5*(d_z+fp1->k_z[n]*alpha_z) ) /
			( fp1->k_z[n] +      theta*dt*0.5*(d_z+fp1->k_z[n]*alpha_z) );
			fp23->bh_z[n] = ( fp1->kh_z[n] - (1.-theta)*dt*0.5*(dh_z+fp1->kh_z[n]*alphah_z) ) /
			( fp1->kh_z[n] +      theta*dt*0.5*(dh_z+fp1->kh_z[n]*alphah_z) );
			fp23->c_z[n] = -d_z * dt*0.5 /
			( fp1->k_z[n]*(fp1->k_z[n]+theta*dt*0.5*(d_z+fp1->k_z[n]*alpha_z)));
			fp23->ch_z[n] = -dh_z * dt*0.5 /
			( fp1->kh_z[n]*(fp1->kh_z[n]+theta*dt*0.5*(dh_z+fp1->kh_z[n]*alphah_z)));
			
			fp4->b_z[n] = ( fp1->k_z[n] - (1.-theta)*dt*(d_z+fp1->k_z[n]*alpha_z) ) /
			( fp1->k_z[n] +      theta*dt*(d_z+fp1->k_z[n]*alpha_z) );
			fp4->bh_z[n] = ( fp1->kh_z[n] - (1.-theta)*dt*(dh_z+fp1->kh_z[n]*alphah_z) ) /
			( fp1->kh_z[n] +      theta*dt*(dh_z+fp1->kh_z[n]*alphah_z) );
			fp4->c_z[n] = -d_z * dt /
			( fp1->k_z[n]*(fp1->k_z[n]+theta*dt*(d_z+fp1->k_z[n]*alpha_z)));
			fp4->ch_z[n] = -dh_z * dt /
			( fp1->kh_z[n]*(fp1->kh_z[n]+theta*dt*(dh_z+fp1->kh_z[n]*alphah_z)));
			
			//printf("\n%zd %lg  %lg  %lg  %lg", n, fp23->c_z[n], fp23->ch_z[n], fp23->b_z[n], fp23->bh_z[n]);
			//printf("\n                                                                 %zd %lg  %lg  %lg  %lg", n, fp4->c_z[n], fp4->ch_z[n], fp4->b_z[n], fp4->bh_z[n]);
		}
		
		for ( size_t n=0; n<g->ab.np+1; ++n ) {
			
			fp1->kH_x[n] = 1.+(g->ab.kappa_max-1.) * pow(((g->ab.np-n+0.5)*g->dx/Lx), g->ab.pmlOrder);
			fp1->kH_z[n] = fp1->kH_x[n];
			
			double dH_x = d0x * pow(((g->ab.np-n+0.5)*g->dx/Lx), g->ab.pmlOrder);
			double alphaH_x = (n-0.5)*da;
			
			fp1->bH_x[n] = 1.0;
			fp1->cH_x[n] = 0.0;
			
			fp23->bH_x[n] = ( fp1->kH_x[n] - (1.-theta)*dt*0.5*(dH_x+fp1->kH_x[n]*alphaH_x) ) /
			( fp1->kH_x[n] +      theta*dt*0.5*(dH_x+fp1->kH_x[n]*alphaH_x) );
			fp23->cH_x[n] = -dH_x * dt*0.5 /
			( fp1->kH_x[n]*(fp1->kH_x[n]+theta*dt*0.5*(dH_x+fp1->kH_x[n]*alphaH_x)));
			
			fp4->bH_x[n] = ( fp1->kH_x[n] - (1.-theta)*dt*(dH_x+fp1->kH_x[n]*alphaH_x) ) /
			( fp1->kH_x[n] +      theta*dt*(dH_x+fp1->kH_x[n]*alphaH_x) );
			fp4->cH_x[n] = -dH_x * dt /
			( fp1->kH_x[n]*(fp1->kH_x[n]+theta*dt*(dH_x+fp1->kH_x[n]*alphaH_x)));
			
			double dH_z = d0z * pow(((g->ab.np-n+0.5)*g->dz/Lz), g->ab.pmlOrder);
			double alphaH_z = (n-0.5)*da;
			
			fp1->bH_z[n] = 1.0;
			fp1->cH_z[n] = 0.0;
			
			fp23->bH_z[n] = ( fp1->kH_z[n] - (1.-theta)*dt*0.5*(dH_z+fp1->kH_z[n]*alphaH_z) ) /
			( fp1->kH_z[n] +      theta*dt*0.5*(dH_z+fp1->kH_z[n]*alphaH_z) );
			fp23->cH_z[n] = -dH_z * dt*0.5 /
			( fp1->kH_z[n]*(fp1->kH_z[n]+theta*dt*0.5*(dH_z+fp1->kH_z[n]*alphaH_z)));
			
			fp4->bH_z[n] = ( fp1->kH_z[n] - (1.-theta)*dt*(dH_z+fp1->kH_z[n]*alphaH_z) ) /
			( fp1->kH_z[n] +      theta*dt*(dH_z+fp1->kH_z[n]*alphaH_z) );
			fp4->cH_z[n] = -dH_z * dt /
			( fp1->kH_z[n]*(fp1->kH_z[n]+theta*dt*(dH_z+fp1->kH_z[n]*alphaH_z)));
			
			//printf("\n%zd %lg  %lg", n, fp23->cH_z[n], fp23->bH_z[n]);
			//printf("\n                                                                 %zd %lg  %lg", n, fp4->cH_z[n], fp4->bH_z[n]);
		}			
	} else {
		for ( size_t n=0; n<g->ab.np; ++n ) {
			
			fp1->k_x[n] = fp1->k_z[n] = fp1->kh_x[n] = fp1->kh_z[n] = 1.;
			
			fp1->b_x[n] = fp1->c_x[n] = fp1->bh_x[n] = fp1->ch_x[n] = 0.;
			fp23->b_x[n] = fp23->bh_x[n] = fp23->c_x[n] = fp23->ch_x[n] = 0.;
			fp4->b_x[n] = fp4->bh_x[n] = fp4->c_x[n] = fp4->ch_x[n] = 0.;
			
			fp1->b_z[n] = fp1->c_z[n] = fp1->bh_z[n] = fp1->ch_z[n] = 0.;
			fp23->b_z[n] = fp23->bh_z[n] = fp23->c_z[n] = fp23->ch_z[n] = 0.;
			fp4->b_z[n] = fp4->bh_z[n] = fp4->c_z[n] = fp4->ch_z[n] = 0.;
		}
		for ( size_t n=0; n<g->ab.np+1; ++n ) {
			
			fp1->kH_x[n] = fp1->kH_z[n] = 1.;
			
			fp1->bH_x[n] = fp1->cH_x[n] = 0.;
			fp23->bH_x[n] = fp23->cH_x[n] = 0.;
			fp4->bH_x[n] = fp4->cH_x[n] = 0.;
			
			fp1->bH_z[n] = fp1->cH_z[n] = 0.;
			fp23->bH_z[n] = fp23->cH_z[n] = 0.;
			fp4->bH_z[n] = fp4->cH_z[n] = 0.;
		}
	}
}


/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */


#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }
double quick_select(double arr[], size_t n) {
	size_t low, high;
	size_t median;
	size_t middle, ll, hh;
	
	low = 0 ; high = n-1 ;
	median = (low + high) / 2;
	for (;;) {
		if (high <= low) /* One element only */
			return arr[median];
		if (high == low + 1) { /* Two elements only */
			if (arr[low] > arr[high])
				ELEM_SWAP(arr[low], arr[high]);
			return arr[median];
		}
	
		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		if (arr[middle] > arr[high]) ELEM_SWAP(arr[middle], arr[high]);
		if (arr[low] > arr[high])    ELEM_SWAP(arr[low], arr[high]);
		if (arr[middle] > arr[low])  ELEM_SWAP(arr[middle], arr[low]);
		
		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP(arr[middle], arr[low+1]);
		
		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
			do ll++; while (arr[low] > arr[ll]);
			do hh--; while (arr[hh] > arr[low]);
			
			if (hh < ll)
				break;
			
			ELEM_SWAP(arr[ll], arr[hh]) ;
		}
		
		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP(arr[low], arr[hh]);
		
		/* Re-set active partition */
		if (hh <= median)
			low = ll;
		
		if (hh >= median)
			high = hh - 1;
	}
}
#undef ELEM_SWAP


void pml_median_average( double *data, const struct grid *g ) {
	
	double *tmp;
	double array[441];
	size_t nmed=1;
	
	size_t nnodes = g->nx2 * g->nz2;
    
    if ( NULL == ( tmp = (double *) malloc( nnodes * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	for ( size_t n=0; n<nnodes; ++n )
		tmp[n] = data[n];
	
	
	for ( ptrdiff_t n=1; n<=g->ab.np; ++n ) {
		
		for ( ptrdiff_t i=(g->ab.np-n); i<g->nx2-(g->ab.np-n); i++ ) {
			ptrdiff_t j=g->ab.np-n;
			
			size_t nel=0;
			for ( ptrdiff_t in=i-nmed; in<=i+nmed; in++ ) {
				for ( ptrdiff_t jn=j-nmed; jn<=j+nmed; jn++ ) {
					if ( in>=0 && in<g->nx2 && jn>=0 && jn<g->nz2 ) {
						array[nel] = data[ in*g->nz2+jn ];						
						nel++;
					}
				}
			}
			tmp[ i*g->nz2+j ] = quick_select(array, nel);
			
			
			j = g->nz2-(g->ab.np-n)-1;
			nel=0;
			for ( ptrdiff_t in=i-nmed; in<=i+nmed; in++ ) {
				for ( ptrdiff_t jn=j-nmed; jn<=j+nmed; jn++ ) {
					if ( in>=0 && in<g->nx2 && jn>=0 && jn<g->nz2 ) {
						array[nel] = data[ in*g->nz2+jn ];						
						nel++;
					}
				}
			}
			tmp[ i*g->nz2+j ] = quick_select(array, nel);
		}
		
		for ( ptrdiff_t j=g->ab.np-n+1; j<g->nz2-g->ab.np-n-1; ++j ) {
			ptrdiff_t i=(g->ab.np-n);
			
			size_t nel=0;
			for ( ptrdiff_t in=i-nmed; in<=i+nmed; in++ ) {
				for ( ptrdiff_t jn=j-nmed; jn<=j+nmed; jn++ ) {
					if ( in>=0 && in<g->nx2 && jn>=0 && jn<g->nz2 ) {
						array[nel] = data[ in*g->nz2+jn ];						
						nel++;
					}
				}
			}
			tmp[ i*g->nz2+j ] = quick_select(array, nel);
		
			i = g->nx2-(g->ab.np-n)-1;
			nel=0;
			for ( ptrdiff_t in=i-nmed; in<=i+nmed; in++ ) {
				for ( ptrdiff_t jn=j-nmed; jn<=j+nmed; jn++ ) {
					if ( in>=0 && in<g->nx2 && jn>=0 && jn<g->nz2 ) {
						array[nel] = data[ in*g->nz2+jn ];						
						nel++;
					}
				}
			}
			tmp[ i*g->nz2+j ] = quick_select(array, nel);
		}		
		
		nmed = (nmed<11) ? nmed++ : nmed;
		for ( size_t n=0; n<nnodes; ++n )
			data[n] = tmp[n];
	}
	
	free( tmp );
}
