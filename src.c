/*
 *  src.c
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
#include <stdio.h>
#include <stdlib.h>

#include "src.h"


void compute_src_fct(struct sourceParams *s, const double dt0, const double fac) {
	
	const double pi = 4.0*atan(1.0);
	const double pi2 = pi*pi;
	double dt = fac*dt0;  // we might need src fct at half time step...
	
	for ( size_t n=0; n<s->nsrc*s->nTemplate; ++n ) {
		double tcut = 1.0/s->s[n].f;
		s->s[n].length = ceil( tcut*2.0/dt );
		if ( NULL == ( s->s[n].fct = malloc( s->s[n].length * sizeof(double))))
		{ fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        double Amax = 0.0;
		for ( size_t i=0; i<s->s[n].length; ++i ) {
			double time = i*dt-tcut;
			double K = pi2 * s->s[n].f * s->s[n].f *time*time;
			s->s[n].fct[i] = s->s[n].A * 2.*pi2*s->s[n].f*s->s[n].f*(1.0-2.0*K)*exp(-K);
            Amax = Amax > fabs(s->s[0].fct[i]) ? Amax : fabs(s->s[0].fct[i]);
		}
        for ( size_t i=0; i<s->s[0].length; ++i ) {
            s->s[0].fct[i] *= s->s[n].A/Amax;
        }
	}
}

void normalize_src_fct(struct sourceParams *s, const struct grid *g) {
	
	double i_area = 1./(g->dx*g->dz);
	for ( size_t n=0; n<s->nsrc*s->nTemplate; ++n ) {
		for ( size_t i=0; i<s->s[n].length; ++i ) {
			s->s[n].fct[i] *= i_area;
		}
	}
}

void add_src_rk4(struct sourceParams *s, double *W, double *Wtmp, const double dt,
				 const double dttmp, const size_t it, const size_t nnodes) {
	
	double *tau_xx    = &(W[0]);
	double *tau_zz    = &(W[nnodes]);
	double *tau_xz    = &(W[2*nnodes]);
	double *p         = &(W[3*nnodes]);
	double *tau_xxtmp = &(Wtmp[0]);
	double *tau_zztmp = &(Wtmp[nnodes]);
	double *tau_xztmp = &(Wtmp[2*nnodes]);
	double *ptmp      = &(Wtmp[3*nnodes]);
	
	for (size_t ns=0; ns<s->nsrc*s->nTemplate; ++ns) {
		s->s[ns].it = 2*it+1;
		if ( (s->s[ns].type == SF || s->s[ns].type == BULK) && s->s[ns].it<s->s[ns].length ) {
			if ( dttmp != 0.0 )
				ptmp[ s->s[ns].i ] += dttmp * s->s[ns].fct[ s->s[ns].it ];
			p[ s->s[ns].i ]    += dt    * s->s[ns].fct[ s->s[ns].it ];
		}
		if ( (s->s[ns].type == SX || s->s[ns].type == BULK || s->s[ns].type == BULK_S) && s->s[ns].it<s->s[ns].length ) {
			if ( dttmp != 0.0 )
				tau_xxtmp[ s->s[ns].i ] += dttmp * s->s[ns].fct[ s->s[ns].it ];
			tau_xx[ s->s[ns].i ]    += dt    * s->s[ns].fct[ s->s[ns].it ];
		}
		if ( (s->s[ns].type == SZ || s->s[ns].type == BULK || s->s[ns].type == BULK_S) && s->s[ns].it<s->s[ns].length ) {
			if ( dttmp != 0.0 )
				tau_zztmp[ s->s[ns].i ] += dttmp * s->s[ns].fct[ s->s[ns].it ];
			tau_zz[ s->s[ns].i ]    += dt    * s->s[ns].fct[ s->s[ns].it ];
		}
		if ( s->s[ns].type == SXZ && s->s[ns].it<s->s[ns].length ) {
			if ( dttmp != 0.0 )
				tau_xztmp[ s->s[ns].i ] += dttmp * s->s[ns].fct[ s->s[ns].it ];
			tau_xz[ s->s[ns].i ]    += dt    * s->s[ns].fct[ s->s[ns].it ];
		}
	}	
}

void add_src_ve_rk4(struct sourceParams *s, double *W, double *Wtmp, const double dt,
                    const double dttmp, const size_t it, const size_t nnodes) {
	
	double *tau_xx    = &(W[0]);
	double *tau_zz    = &(W[nnodes]);
	double *tau_xz    = &(W[2*nnodes]);
	double *tau_xxtmp = &(Wtmp[0]);
	double *tau_zztmp = &(Wtmp[nnodes]);
	double *tau_xztmp = &(Wtmp[2*nnodes]);
	
	for (size_t ns=0; ns<s->nsrc*s->nTemplate; ++ns) {
		s->s[ns].it = 2*it+1;
		if ( (s->s[ns].type == SX || s->s[ns].type == BULK) && s->s[ns].it<s->s[ns].length ) {
			if ( dttmp != 0.0 )
				tau_xxtmp[ s->s[ns].i ] += dttmp * s->s[ns].fct[ s->s[ns].it ];
			tau_xx[ s->s[ns].i ]    += dt    * s->s[ns].fct[ s->s[ns].it ];
		}
		if ( (s->s[ns].type == SZ || s->s[ns].type == BULK) && s->s[ns].it<s->s[ns].length ) {
			if ( dttmp != 0.0 )
				tau_zztmp[ s->s[ns].i ] += dttmp * s->s[ns].fct[ s->s[ns].it ];
			tau_zz[ s->s[ns].i ]    += dt    * s->s[ns].fct[ s->s[ns].it ];
		}
		if ( s->s[ns].type == SXZ && s->s[ns].it<s->s[ns].length ) {
			if ( dttmp != 0.0 )
				tau_xztmp[ s->s[ns].i ] += dttmp * s->s[ns].fct[ s->s[ns].it ];
			tau_xz[ s->s[ns].i ]    += dt    * s->s[ns].fct[ s->s[ns].it ];
		}
	}
}

void add_src_ve_sh_rk4(struct sourceParams *s, double *W, double *Wtmp, const double dt,
					const double dttmp, const size_t it, const size_t nnodes) {
	
	double *tau_xy    = &(W[0]);
	double *tau_yz    = &(W[nnodes]);
	double *tau_xytmp = &(Wtmp[0]);
	double *tau_yztmp = &(Wtmp[nnodes]);
	
	for (size_t ns=0; ns<s->nsrc*s->nTemplate; ++ns) {
		s->s[ns].it = 2*it+1;
		
		// TODO: check this
		
		if ( (s->s[ns].type == SXY) && s->s[ns].it<s->s[ns].length ) {
			if ( dttmp != 0.0 )
				tau_xytmp[ s->s[ns].i ] += dttmp * s->s[ns].fct[ s->s[ns].it ];
			tau_xy[ s->s[ns].i ]    += dt    * s->s[ns].fct[ s->s[ns].it ];
		}
		if ( (s->s[ns].type == SYZ) && s->s[ns].it<s->s[ns].length ) {
			if ( dttmp != 0.0 )
				tau_yztmp[ s->s[ns].i ] += dttmp * s->s[ns].fct[ s->s[ns].it ];
			tau_yz[ s->s[ns].i ]    += dt    * s->s[ns].fct[ s->s[ns].it ];
		}
	}
}

void add_src(const struct sourceParams *s, double *txx, double *tzz, double *txz,
			 const double dt, const size_t it) {
    
	for (size_t ns=0; ns<s->nsrc*s->nTemplate; ++ns) {
		if ( (s->s[ns].type == SX || s->s[ns].type == BULK) && it < s->s[ns].length ) {
			txx[ s->s[ns].i ]    += dt * s->s[ns].fct[ it ];
		}
		if ( (s->s[ns].type == SZ || s->s[ns].type == BULK) && it < s->s[ns].length ) {
			tzz[ s->s[ns].i ]    += dt * s->s[ns].fct[ it ];
		}
		if ( s->s[ns].type == SXZ && it < s->s[ns].length ) {
			txz[ s->s[ns].i ]    += dt * s->s[ns].fct[ it ];
		}
	}
}

void add_force_src_cyl(const struct sourceParams *s, struct variables_cyl *v,
                       const double dt, const size_t it) {
    
    for (size_t ns=0; ns<s->nsrc*s->nTemplate; ++ns) {
        if ( (s->s[ns].type == FR || s->s[ns].type == BULK) && it < s->s[ns].length ) {
            v->vr[ s->s[ns].i ] += dt*v->bj[ s->s[ns].i ] * s->s[ns].fct[ it ];
        }
        if ( s->s[ns].type == FT && it < s->s[ns].length ) {
            v->vt[ s->s[ns].i ] += dt*v->bij[ s->s[ns].i ] * s->s[ns].fct[ it ];
        }
        if ( (s->s[ns].type == FZ || s->s[ns].type == BULK) && it < s->s[ns].length ) {
            v->vz[ s->s[ns].i ] += dt*v->bi[ s->s[ns].i ] * s->s[ns].fct[ it ];
        }
    }
}

void compute_kurkjian(struct sourceParams *s, const struct variables_cyl *v,
                      const struct grid *g, const double dt, const int n) {
    
    const double pi = 4.0*atan(1.0);
    const double pi2 = pi*pi;

    double epsilon_m = 1.0;
    if ( n>0 ) epsilon_m = 2.0;
    double term_kurkjian = v->lij[s->s[0].i]/(2.*pi)*epsilon_m *
    2./(g->dx*g->dx*g->dz) * pow(s->s[0].x/g->dx, n);
    
    
    double tcut = 1.0/s->s[0].f;
    s->s[0].length = ceil( tcut*2.0/dt );
    if ( NULL == ( s->s[0].fct = malloc( s->s[0].length * sizeof(double)))) {
        fprintf(stderr, "Error: cannot allocate memory\n"); abort();
    }
    double Amax = 0.0;
    for ( size_t i=0; i<s->s[0].length; ++i ) {
        double time = i*dt-tcut;
        double K = pi2 * s->s[0].f * s->s[0].f *time*time;
        s->s[0].fct[i] = 2.*pi2*s->s[0].f*s->s[0].f*(1.0-2.0*K)*exp(-K);
        Amax = Amax > fabs(s->s[0].fct[i]) ? Amax : fabs(s->s[0].fct[i]);
    }
    for ( size_t i=0; i<s->s[0].length; ++i ) {
        s->s[0].fct[i] *= s->s[0].A/Amax * term_kurkjian;
    }
}

void add_kurkjian_src(const struct sourceParams *s, struct variables_cyl *v,
                      const double dt, const size_t it) {
    if ( it < s->s[0].length ) {
        v->trr[ s->s[0].i ] += dt*s->s[0].fct[ it ];
        v->tzz[ s->s[0].i ] += dt*s->s[0].fct[ it ];
    }
}
