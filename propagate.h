/*
 *  propagate.h
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


#ifndef __propagate_H__
#define __propagate_H__

#include "structs.h"

void init_fftw_data(double *, struct fftw_data *, const struct grid *);
void free_fftw_data(struct fftw_data *);

void init_fftw_data_ve(double *, struct fftw_data *, const struct grid *);
void free_fftw_data_ve(struct fftw_data *);
void init_fftw_data_ve_sh(double *, struct fftw_data_ve_sh *, const struct grid *);
void free_fftw_data_ve_sh(struct fftw_data_ve_sh *);

void propagate(double *, double *, const struct grid *,
               const struct computationVariables *, struct fftw_data *);

void propagateCPML(double *, double *, const struct grid *,
                   const struct computationVariables *, struct fftw_data *,
                   struct mem_pml *, struct fac_pml *);

void propagateVTI(double *, double *, const struct grid *,
				  const struct computationVariablesVTI *, struct fftw_data *);

void propagateVTI_CPML(double *, double *, const struct grid *,
					   const struct computationVariablesVTI *,
					   struct fftw_data *, struct mem_pml *,
					   struct fac_pml *);

void propagateVE_CPML(double *, double *, const struct grid *,
					   const struct computationVariablesVE_VTI *,
					   struct fftw_data *, struct mem_pml *,
					   struct fac_pml *, const int);

void propagateVE_SH_CPML(double *, double *, const struct grid *,
						 const struct computationVariablesVE_SH_VTI *,
						 struct fftw_data_ve_sh *, struct mem_pml_sh *,
						 struct fac_pml *);

void partial_x_forward(struct fftw_data *, const struct grid *);
void partial_x_backward(struct fftw_data *, const struct grid *);
void partial_z_forward(struct fftw_data *, const struct grid *);
void partial_z_backward(struct fftw_data *, const struct grid *);

void partial_x_forward_sh(struct fftw_data_ve_sh *, const struct grid *);
void partial_x_backward_sh(struct fftw_data_ve_sh *, const struct grid *);
void partial_z_forward_sh(struct fftw_data_ve_sh *, const struct grid *);
void partial_z_backward_sh(struct fftw_data_ve_sh *, const struct grid *);

void update_vx(double *, const double *, const double *, const double *,
               const double *, const double *, const double *, const double *,
               const struct grid *, const double, const double,
               const size_t);
void update_vz(double *, const double *, const double *, const double *,
               const double *, const double *, const double *, const double *,
               const struct grid *, const double, const double,
               const size_t);

void update_txxzz(double *, double *, const double *, const double *,
                  const double *, const double *,  const double *, const double *, 
                  const double *, const double *,
                  const struct grid *, const double, const double,
                  const size_t);
void update_txz(double *, const double *, const double *, const double *, 
                const double *, const double *, const double *, const double *,
                const struct grid *, const double, const double,
                const size_t);

void compute_div(double *, const struct grid *, struct fftw_data *);
void compute_curl(double *, const struct grid *, struct fftw_data *);
#endif
