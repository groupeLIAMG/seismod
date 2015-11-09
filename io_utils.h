/*
 *  io_utils.h
 *
 *  Created by Bernard Giroux on 10-08-28.
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

#ifndef __IO_UTILS_H__
#define __IO_UTILS_H__

#include "structs.h"

void set_defaults(struct grid *g, struct inputParams *p);
void set_defaults_3d(struct grid3d *g, struct inputParams *p);
void print_usage(char *, FILE *out, int exit_code);
void process_args(int argc, char * const argv[], struct inputParams *p);
void read_grid_params(const char filename[], struct grid *g);
void read_grid_params_ve(const char filename[], struct grid *g,
                         struct materialPropertiesVE_VTI *mp);
void read_grid_params_ve_sh(const char filename[], struct grid *g,
							struct materialPropertiesVE_SH_VTI *mp);
void read_grid3d_params(const char filename[], struct grid3d *g);
void read_model(const char filename[], const struct grid *g,
                struct materialProperties *mp);
void read_modelVTI(const char filename[], const struct grid *g, 
				   struct materialPropertiesVTI *mp);
void read_model_e(const char filename[], const struct grid *g, double *vp,
				  double *vs, double *rho);
void read_model_e_cyl(const char filename[], const struct grid *g, double *vp,
                      double *vs, double *rho);
void read_model_a(const char filename[], const struct grid3d *g, double *vp,
				  double *rho);
void read_model_ve(const char filename[], const struct grid *g,
				   struct materialPropertiesVE_VTI *mp);
void read_model_ve_sh(const char filename[], const struct grid *g,
					  struct materialPropertiesVE_SH_VTI *mp);
void read_source(const char filename[], const struct grid *, struct sourceParams *);
void read_source_3d(const char filename[], const struct grid3d *, struct sourceParams *);
void read_output(const char filename[], const struct grid *, struct outputParams *);
void read_output_3d(const char filename[], const struct grid3d *, struct outputParams *);
void write_trace(const double *, const double, struct outputParams *, const size_t);
void fill_data_segy(const double *, const double, struct outputParams *, const size_t);
void write_snapshot(const double *, const double, const struct grid *,
                    struct outputParams *, const size_t, const short);
void write_snapshot_nc(const double *, const double, const struct grid *,
                       struct outputParams *, const size_t, const short);
void write_snapshot3d_nc(const double *, const double, const struct grid3d *,
                       struct outputParams *, const size_t, const short);
void write_field_nc(const double *, const char *, const char *,
                    const struct grid *, struct outputParams *, const short);
void save_checkpoint(const size_t, const struct grid *,
					 const struct inputParams *, const struct sourceParams *,
					 const struct outputParams *, const struct saveEnergy *,
					 const struct fac_pml *, const struct fac_pml *,
					 const struct fac_pml *, const struct mem_pml mem[3],
					 const struct computationVariables *,
					 const double *, const double *, const double *,
					 const double *, const double *, const double *,
					 const double *, const double *	);
void read_checkpoint1(const char *, size_t *, struct grid *,
					  struct inputParams *);
void read_checkpoint2(const char *, const struct grid *,
					  const struct inputParams *, struct sourceParams *,
					  struct outputParams *, struct saveEnergy *,
					  struct fac_pml *, struct fac_pml *, struct fac_pml *,
					  struct mem_pml mem[3], struct computationVariables *,
					  double *, double *, double *, double *, double *,
					  double *, double *, double *);
void save_checkpointVTI(const size_t, const struct grid *,
						const struct inputParams *, const struct sourceParams *,
						const struct outputParams *, const struct saveEnergyVTI *,
						const struct fac_pml *, const struct fac_pml *,
						const struct fac_pml *, const struct mem_pml mem[3],
						const struct computationVariablesVTI *,
						const double *, const double *, const double *,
						const double *, const double *, const double *,
						const double *, const double *	);
void read_checkpoint2VTI(const char *, const struct grid *,
						 const struct inputParams *, struct sourceParams *,
						 struct outputParams *, struct saveEnergyVTI *,
						 struct fac_pml *, struct fac_pml *, struct fac_pml *,
						 struct mem_pml mem[3], struct computationVariablesVTI *,
						 double *, double *, double *, double *, double *,
						 double *, double *, double *);
void save_segy(const struct inputParams *, const struct outputParams *,
               const struct sourceParams *);

void check_dt_trc(struct outputParams *, double);

#endif
