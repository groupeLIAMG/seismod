/*
 *  src.h
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

#ifndef __src_H__
#define __src_H__

#include "structs.h"

void compute_src_fct(struct sourceParams *, const double, const double);
void normalize_src_fct(struct sourceParams *, const struct grid *);
void add_src_rk4(struct sourceParams *, double *, double *,
				 const double, const double, const size_t, const size_t);
void add_src_ve_rk4(struct sourceParams *, double *, double *,
                    const double, const double, const size_t, const size_t);
void add_src_ve_sh_rk4(struct sourceParams *, double *, double *,
                    const double, const double, const size_t, const size_t);
void add_src(const struct sourceParams *, double *, double *, double *,
			 const double, const size_t);
void add_force_src_cyl(const struct sourceParams *s, struct variables_cyl *v,
                       const double dt, const size_t it);
void compute_kurkjian(struct sourceParams *s, const struct variables_cyl *v,
                      const struct grid *g, const double dt, const int n);
void add_kurkjian_src(const struct sourceParams *s, struct variables_cyl *v,
                      const double dt, const size_t it);

#endif
