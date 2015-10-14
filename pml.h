/*
 *  pml.h
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


#include <stdlib.h>

#include "structs.h"

void pml_v(double *, double *, double *, double *, const double *, const double *,
		   const double *, const double *, const double *, const double *,
		   const double [], const double [], const double [],
		   const double *, const double *, const double *,
		   const double *, const double *,
		   const size_t, const size_t, const size_t, const double, const double);

void pml_t(const double *, const double *, const double *, const double *,
           double *, double *, double *, double *, double *, double *,
		   const double [], const double [], const double [],
		   const double *, const double *,
		   const double *, const double *, const double *,
		   const size_t, const size_t, const size_t, const double, const double);

void alloc_cpml(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
				struct mem_pml mem[3], const struct grid *g);

void free_cpml(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
			   struct mem_pml mem[3]);

void alloc_cpml_ve(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
                   struct mem_pml mem[3], const struct grid *g);

void free_cpml_ve(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
                  struct mem_pml mem[3]);

void alloc_cpml_ve_sh(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
					  struct mem_pml_sh mem_sh[3], const struct grid *g);

void free_cpml_ve_sh(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
					 struct mem_pml_sh mem[3]);

void compute_cpml(struct fac_pml *fp1, struct fac_pml *fp23, struct fac_pml *fp4,
				  const struct grid *g, const double dt, const double f,
				  const short iwipe);
double quick_select(double arr[], size_t n);
void pml_median_average( double *data, const struct grid *g );