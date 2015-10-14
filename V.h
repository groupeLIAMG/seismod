/*
 * Copyright (c) 2012, Bernard Giroux
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


#ifndef _V_H_
#define _V_H_

#include <complex.h>


double complex VqP1(const double lx,
                    const double lz,
                    const double cu11,
                    const double cu13,
                    const double cu33,
                    const double cu55,
                    const double alpha1,
                    const double alpha3,
                    const double M,
                    const double rho,
                    const double rho_f,
                    const double nk1,
                    const double nk3,
                    const double m1,
                    const double m3,
                    const double omega);

double complex VqP2(const double lx,
                    const double lz,
                    const double cu11,
                    const double cu13,
                    const double cu33,
                    const double cu55,
                    const double alpha1,
                    const double alpha3,
                    const double M,
                    const double rho,
                    const double rho_f,
                    const double nk1,
                    const double nk3,
                    const double m1,
                    const double m3,
                    const double omega);


double complex VqS(const double lx,
                   const double lz,
                   const double cu11,
                   const double cu13,
                   const double cu33,
                   const double cu55,
                   const double alpha1,
                   const double alpha3,
                   const double M,
                   const double rho,
                   const double rho_f,
                   const double nk1,
                   const double nk3,
                   const double m1,
                   const double m3,
                   const double omega);


#endif
