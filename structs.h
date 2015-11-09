/*
 *  structs.h
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


#ifndef __STRUCTS_H__
#define __STRUCTS_H__

#include <complex.h>
#include <stddef.h>

#include "fftw3.h"

enum sourceType { SX, SY, SZ, SXY, SXZ, SYZ, BULK, SF, BULK_S, FR, FT, FZ, KURKJIAN };
enum component { VX, VY, VZ, QX, QZ, TXX, TZZ, TXZ, P, TXY, TYZ, DIV, CURL, VR, VT, TRR, TRT, TRZ, TTT, TTZ };
enum typeRecord { TRACE, SNAPSHOT };
enum coordinates { CARTESIAN, CYLINDRICAL };

struct abs_params {
    size_t np;           // number of padding nodes for absorbing boundaries
	short  median;       // smooth absorbing region with median filter
	short  alpha_max;    // if 0, alpha = 0.0; if 1, alpha_max = pi * f
    double alfa;
    double vmax;
	double pmlOrder;
	double Rc;           // PML, theoretical coefficient of reflection
	double kappa_max;    // kappa CPML (Komatitsch et Martin, 2007)
    double *gobx;
    double *gobz;
};

struct grid {
    enum coordinates coord;
    size_t nx, nz;     // number of nodes in x and z
    size_t nx2, nz2;   // total number of nodes
    short  Qdamping;   // decrease Q in absorbing boundaries (1: yes, 0: no)
    double x0, z0;     // origin                                           [ m ]
    double dx, dz;     // grid cell size in x and z                        [ m ]
    struct abs_params ab;
};

struct grid3d {
    size_t nx, ny, nz;      // number of nodes in x and z
    size_t nx2, ny2, nz2;   // total number of nodes
    double x0, y0, z0;      // origin                                      [ m ]
    double dx, dy, dz;      // grid cell size in x and z                   [ m ]
    struct abs_params ab;
};

struct inputParams {
    double dt;         //                                                 [ ms ]
    double duration;   //                                                 [ ms ]
    char   modelfile[80];
    char   sourcefile[80];
    char   outputfile[80];
    char   basename[80];
	char   checkpointfile[80];
    short  iwipe;
    short  plotStrips;
	short  checkpoint;
    short  segy;
    short  check_model;
    short  hasDiv;
    short  hasCurl;
	int    chkpt_inc;
    struct abs_params *ab;
	short  saveEnergy;
	double roiExmin;
	double roiExmax;
	double roiEzmin;
	double roiEzmax;
    int    shotpt_no;  // shotpoint number (for SEGY files)
    int    simulateCMP;
    int    n;
};

struct materialProperties {
    
    //
    //  Properties of materials
    //
    //  the following are input data, with corresponding units in brakets
    //     note: 1 cP = 1e-3 Pa.s ; 1 mD = 1e-15 m^2
    
    double *K_m;       // bulk modulus of drained matrix                 [ GPa ]
    double *K_s;       // bulk modulus of the solid                      [ GPa ]
    double *K_f;       // bulk modulus of the fluid                      [ GPa ]
    double *phi;       // porosity                                         [ - ]
    double *mu;        // shear modulus of the matrix                    [ GPa ]
    double *rho_s;     // solid density                               [ kg/m^3 ]
    double *rho_f;     // fluid density                               [ kg/m^3 ]
    double *T;         // tortuosity                                       [ - ]
    double *eta;       // fluid viscosity                                 [ cP ]
    double *kappa;     // permeability                                    [ mD ]
    double *Q;         // seismic quality factor                           [ - ]
	double *f0;        // relaxation frequency                            [ Hz ]
};    

struct materialPropertiesVTI {
    
    //
    //  Properties of materials
    //
    //  the following are input data, with corresponding units in brakets
    //     note: 1 cP = 1e-3 Pa.s ; 1 mD = 1e-15 m^2
    
    double *K_m;       // bulk modulus of drained matrix                 [ GPa ]
    double *K_s;       // bulk modulus of the solid                      [ GPa ]
    double *K_f;       // bulk modulus of the fluid                      [ GPa ]
    double *phi;       // porosity                                         [ - ]
    double *mu;        // shear modulus of the matrix                    [ GPa ]
    double *rho_s;     // solid density                               [ kg/m^3 ]
    double *rho_f;     // fluid density                               [ kg/m^3 ]
	double *T1;        // tortuosity (along x_1)                           [ - ]
    double *T3;        // tortuosity (along x_3)                           [ - ]
    double *eta;       // fluid viscosity                                 [ cP ]
	double *kappa1;    // permeability (along x_1)                        [ mD ]
    double *kappa3;    // permeability (along x_3)                        [ mD ]
    double *Q ;        // seismic quality factor                           [ - ]
	double *f0;        // relaxation frequency                            [ Hz ]
	double *epsilon;   // Thomsen anisotropy parameter                     [ - ]
	double *delta;     // Thomsen anisotropy parameter                     [ - ]
//	double *gamma;     // Thomsen anisotropy parameter                     [ - ]
};    

struct materialPropertiesVE_VTI {
    int L;             // Nb of relaxation mechanisms, q-dilatational mode [ - ]
	double *K;         // bulk modulus                                   [ GPa ]
	double *mu;        // shear modulus                                  [ GPa ]
    double *rho;       // density                                     [ kg/m^3 ]
	double *epsilon;   // Thomsen anisotropy parameter                     [ - ]
	double *delta;     // Thomsen anisotropy parameter                     [ - ]
	double **Q1;       // Q for quasi-dilatational mode                    [ - ]
    double **f1;       // Relaxation frequencies - Q1                     [ Hz ]
    double *Q2;        // Q for shear mode                                 [ - ]
    double *f2;        // Relaxation frequency - Q2                       [ Hz ]
};

struct materialPropertiesVE_SH_VTI {
	double *mu;        // shear modulus                                  [ GPa ]
    double *rho;       // density                                     [ kg/m^3 ]
	double *gamma;     // Thomsen anisotropy parameter                     [ - ]
    double *Q;         // Q for shear mode                                 [ - ]
    double *f;         // Relaxation frequency - Q                        [ Hz ]
};


struct computationVariables {
  
	double dt;         // time step
    double *mu;        // shear modulus of the matrix
    
    //
    // intermediate variables, equation numbers refer to Carcione and Helle (1999)
    //
    double *epsilon;   // eq (5)
    double *E;         // stiffness, P-wave modulus of dry skeleton, eq (6)
    double *M;         // coupling modulus, eq (7)
    double *alpha;     // poroelastic coefficient of effective stress, eq (9)
    double *varphi;    // eq (19)
    double *tau_s;     // relaxation time
    
    //
    // variables interpolated on the staggered grid
    //
    double *rho_i;     // rho interpolated at (i+1/2, j)
    double *rho_j;     // rho interpolated at (i, j+1/2)
    double *rho_f_i;   // rho_f interpolated at (i+1/2, j)
    double *rho_f_j;   // rho_f interpolated at (i, j+1/2)
    double *nk_i;      // eta/kappa interpolated at (i+1/2, j)
    double *nk_j;      // eta/kappa interpolated at (i, j+1/2)
    double *mu_ij;     // mu interpolated at (i+1/2, j+1/2)
    double *m_i;       // m interpolated at (i+1/2, j)
    double *m_j;       // m interpolated at (i, j+1/2)
    
};

struct computationVariablesVTI {
	
	double dt;         // time step
    
	double *c11;       // at (i, j)
	double *c13;       // at (i, j)
	double *c33;       // at (i, j)
	double *c55;       // at (i+1/2, j+1/2)
	
    //
    // intermediate variables
    //
    double *epsilon;
    double *M;          // coupling modulus
    double *alpha1;     // poroelastic coefficient of effective stress
    double *alpha3;     // poroelastic coefficient of effective stress
    double *varphi;     // 
    double *tau_s;      // relaxation time
    
    //
    // variables interpolated on the staggered grid
    //
    double *rho_i;     // rho interpolated at (i+1/2, j)
    double *rho_j;     // rho interpolated at (i, j+1/2)
    double *rho_f_i;   // rho_f interpolated at (i+1/2, j)
    double *rho_f_j;   // rho_f interpolated at (i, j+1/2)
	
    double *nk1;       // eta/kappa1 interpolated at (i+1/2, j)
    double *nk3;       // eta/kappa3 interpolated at (i, j+1/2)
    double *m1;        // m1 interpolated at (i+1/2, j)
    double *m3;        // m3 interpolated at (i, j+1/2)
    
};


struct computationVariablesVE_VTI {
	
	double dt;         // time step
    
	double *c11;       // at (i, j)
	double *c13;       // at (i, j)
	double *c33;       // at (i, j)
	double *c55;       // at (i+1/2, j+1/2)
	
    //
    // intermediate variables
    //
    double *K0;
    double *c55_0;
    double *c55_0_ij;
    
    double **epsilon1;
    double **tau_s_1;
    double **eta1;

    double *epsilon2;
    double *epsilon3;
    double *tau_s_2;
    double *tau_s_2_ij;
    
    double *eta2;
    double *eta2_ij;
    
    //
    // variables interpolated on the staggered grid
    //
    double *rho_i;     // rho interpolated at (i+1/2, j)
    double *rho_j;     // rho interpolated at (i, j+1/2)
};

struct computationVariablesVE_SH_VTI {
	
	double dt;         // time step
    
	double *c44;       // at (i, j+1/2)
	double *c66;       // at (i+1/2, j)
	
    //
    // intermediate variables
    //
    double *c44_0;
    double *c66_0;
    
    double *epsilon1;
    double *tau_s_1;
    double *eta1;
	
    double *epsilon2;
    double *tau_s_2;
    double *eta2;
    
	double *rho;
};


struct source {
	enum sourceType type;
	double f;          // frequency of source wavelet                    [ Hz ]
	double A;          // source strength                                [ MPa ]
	double x;          // position along x                               [ m ]
	double y;          // position along y                               [ m ]
	double z;          // position along z                               [ m ]
	size_t i;          // index in grid
	double *fct;       // source function
	size_t it;         // index of time function
	size_t length;     // length of src fct (in terms of sample)
};

struct sourceParams {
	size_t nsrc;
	size_t nTemplate;
	struct source *s;
};

struct fftw_data {
    
    complex double *kx_f, *kx_b, *kz_f, *kz_b;
    fftw_complex *o_x, *o_z;
    double *t1, *t2, *t3;
    fftw_plan tauxx_x_f, tauxx_x_i;
    fftw_plan tauxz_x_f, tauxz_x_i, tauxz_z_f, tauxz_z_i;
    fftw_plan tauzz_z_f, tauzz_z_i;
    fftw_plan p_x_f, p_x_i, p_z_f, p_z_i;
    fftw_plan vx_x_f, vx_x_i, vx_z_f, vx_z_i;
    fftw_plan vz_x_f, vz_x_i, vz_z_f, vz_z_i;
    fftw_plan qx_x_f, qx_x_i;
    fftw_plan qz_z_f, qz_z_i;        
};

struct fftw_data_ve {
    
    complex double *kx_f, *kx_b, *kz_f, *kz_b;
    fftw_complex *o_x, *o_z;
    double *t1, *t2, *t3;
    fftw_plan tauxx_x_f, tauxx_x_i;
    fftw_plan tauxz_x_f, tauxz_x_i, tauxz_z_f, tauxz_z_i;
    fftw_plan tauzz_z_f, tauzz_z_i;
    fftw_plan vx_x_f, vx_x_i, vx_z_f, vx_z_i;
    fftw_plan vz_x_f, vz_x_i, vz_z_f, vz_z_i;
};

struct fftw_data_ve_sh {
    
    complex double *kx_f, *kx_b, *kz_f, *kz_b;
    fftw_complex *o_x, *o_z;
    double *t1, *t2;
    fftw_plan tauxy_x_f, tauxy_x_i;
    fftw_plan tauyz_z_f, tauyz_z_i;
    fftw_plan vy_z_f, vy_z_i;
    fftw_plan vy_x_f, vy_x_i;
};

struct record {
    double x;
    double y;
    double z;
    double dt;
    double t0;
    size_t i;
    enum component comp;
    enum typeRecord type;
    float *data;
    FILE *fid;
};

struct outputParams {
    size_t nrec;
    struct record *r;
    char   basename[80];
};

struct fac_pml {
    double *k_x;  // kappa of Roden and Gedney, 2000
    double *k_z;
    double *kh_x;  // at half grid cell
    double *kh_z;  // at half grid cell
    double *kH_x;  // at half grid cell, end of grid
    double *kH_z;  // at half grid cell, end of grid
    double *b_x;
    double *c_x;
    double *bh_x;  // at half grid cell
    double *ch_x;  // at half grid cell
    double *bH_x;  // at half grid cell, end of grid
    double *cH_x;  // at half grid cell, end of grid
    double *b_z;
    double *c_z;
    double *bh_z;  // at half grid cell
    double *ch_z;  // at half grid cell
    double *bH_z;  // at half grid cell, end of grid
    double *cH_z;  // at half grid cell, end of grid
};

struct fac_cpml_cyl {
    double *ik_r;  // 1/kappa of Roden and Gedney, 2000
    double *ik_z;
    double *ikh_z;  // at half grid cell
    double *ikH_r;  // at half grid cell, end of grid
    double *ikH_z;  // at half grid cell, end of grid
    double *b_r;
    double *c_r;
    double *bH_r;  // at half grid cell, end of grid
    double *cH_r;  // at half grid cell, end of grid
    double *b_z;
    double *c_z;
    double *bh_z;  // at half grid cell
    double *ch_z;  // at half grid cell
    double *bH_z;  // at half grid cell, end of grid
    double *cH_z;  // at half grid cell, end of grid
};

struct mem_pml {  // memory variables
    double *dx_txx;
    double *dx_p;
    double *dx_txz;
    double *dx_vx;
    double *dx_qx;
    double *dx_vz;
    double *dz_txz;
    double *dz_tzz;
    double *dz_p;
    double *dz_vz;
    double *dz_qz;
    double *dz_vx;
};

struct mem_pml_sh {  // memory variables
    double *dx_txy;
    double *dx_vy;
    double *dz_tyz;
    double *dz_vy;
};

struct mem_cpml_cyl {
    // for v_r
    double *dtrr_dr;
    double *trr_r1;
    double *trt_r1;
    double *ttt_r1;
    double *dtrz_dz;
    // for v_t
    double *trt_r2;
    double *dtrt_dr;
    double *ttt_r2;
    double *dttz_dz;
    // for v_z
    double *dtzz_dz;
    double *dtrz_dr;
    double *trz_r;
    double *ttz_r;
    // for \tau_rr, \tau_\theta\theta, \tau_zz
    double *dvr_dr;
    double *vr_r1;
    double *vt_r1;
    double *dvz_dz;
    // for \tau_rz
    double *dvr_dz;
    double *dvz_dr;
    // for \tau_r\theta
    double *dvt_dr;
    double *vr_r2;
    double *vt_r2;
    // for \tau_\theta z
    double *dvt_dz;
    double *vz_r;
};

struct variables_cyl {
    // particle velocities
    double *vr;   // v_r
    double *vt;   // v_\theta
    double *vz;   // v_z
    // stresses
    double *trr;  // \tau_{rr}
    double *ttt;  // \tau_{\theta\theta}
    double *tzz;  // \tau_{zz}
    double *trz;  // \tau_{rz}
    double *trt;  // \tau_{r\theta}
    double *ttz;  // \tau_{\theta z}
    // properties
    double *lij;  // \lambda at (i+1/2,j+1/2)
    double *l2mij;// \lambda+2\mu at (i+1/2,j+1/2)
    double *mu;   // \mu at (i,j)
    double *mui;  // \mu at (i+1/2)
    double *muj;  // \mu at (i,j+1/2)
    double *bi;   // 1/\rho at (i+1/2,j)
    double *bj;   // 1/\rho at (i,j+1/2)
    double *bij;  // 1/\rho at (i+1/2,j+1/2)
};

struct saveEnergy {
    double *rho11;
    double *rho12;
    double *rho22;
	size_t i1E;
	size_t i2E;
	size_t j1E;
	size_t j2E;
	FILE *fid;
};	

struct saveEnergyVTI {
	double *fE1;
    double *fE2;
    double *rEx;
    double *rEz;	
	size_t i1E;
	size_t i2E;
	size_t j1E;
	size_t j2E;
	FILE *fid;
};	


#endif