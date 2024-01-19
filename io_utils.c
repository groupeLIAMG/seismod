/*
 *  io_utils.c
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

#include <sys/stat.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <netcdf.h>
/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}



#include "io_utils.h"

extern int verbose;
extern char *parFormat;
extern char *modelFormat;
extern char *srcFormat;
extern char *rcvFormat;
extern char *optarg;

void set_defaults(struct grid *g, struct inputParams *p) {

    g->coord = CARTESIAN;
    g->nx = 300;
    g->nz = 360;
    g->x0 = 0.0;
    g->z0 = 0.0;
    g->dx = 5.0;
    g->dz = 5.0;
    g->ab.np = 20;
    g->ab.alfa = 0.16;
    g->ab.vmax = 5000.0;
	g->ab.pmlOrder = 2.;
	g->ab.Rc = 1.e-5;
	g->ab.kappa_max = 1.;
	g->ab.alpha_max = 1;
    g->Qdamping = 0;
    
    p->iwipe = 0;
    p->plotStrips = 0;
	p->checkpoint = 0;
    p->segy = 0;
    p->check_model = 0;
	p->chkpt_inc = 0;
    p->dt = 1.e-4;
    p->duration = 1.2;
	strcpy(p->basename, "run01");
	p->saveEnergy = 0;
	p->roiExmin = 0.0;
	p->roiExmax = 0.0;
	p->roiEzmin = 0.0;
	p->roiEzmax = 0.0;
    p->shotpt_no = 1;
    p->simulateCMP = 0;
    p->n = 0;
}

void set_defaults_3d(struct grid3d *g, struct inputParams *p) {
    
    g->nx = 300;
    g->ny = 300;
    g->nx = 360;
    g->x0 = 0.0;
    g->y0 = 0.0;
    g->z0 = 0.0;
    g->dx = 5.0;
    g->dy = 5.0;
    g->dz = 5.0;
    g->ab.np = 20;
    g->ab.alfa = 0.16;
    g->ab.vmax = 5000.0;
	g->ab.pmlOrder = 2.;
	g->ab.Rc = 1.e-5;
	g->ab.kappa_max = 1.;
	g->ab.alpha_max = 1;
    
    p->iwipe = 0;
    p->plotStrips = 0;
	p->checkpoint = 0;
    p->segy = 0;
    p->check_model = 0;
	p->chkpt_inc = 0;
    p->dt = 1.e-4;
    p->duration = 1.2;
	strcpy(p->basename, "run01");
	p->saveEnergy = 0;
	p->roiExmin = 0.0;
	p->roiExmax = 0.0;
	p->roiEzmin = 0.0;
	p->roiEzmax = 0.0;
    p->shotpt_no = 1;
    p->simulateCMP = 0;
}

void print_usage( char *prog, FILE *out, int exit_code ) {
    fprintf(out, "Usage: %s [options] -p parameter_file.dat\n", prog);
	fprintf(out, "\noptions are:\n");
    fprintf(out, "  -h    Print this message\n");
    fprintf(out, "  -v    Verbose mode (type twice for increased verbosity)\n");
	fprintf(out, "  -c checkpoint_file.dat   Restart run using data in checkpoint_file.dat (-p argument ignored)\n\n");
    fprintf(out, "\n /*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\ \n\n");
    fprintf(out, "%s", parFormat);
    fprintf(out, "\n /*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\ \n\n");
    fprintf(out, "%s", modelFormat);
    fprintf(out, "\n /*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\ \n\n");
    fprintf(out, "%s", srcFormat);
    fprintf(out, "\n /*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\ \n\n");
    fprintf(out, "%s", rcvFormat);
    fprintf(out, "\n /*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\*/*\\ \n\n");
	fflush(out);
    exit( exit_code );
}

void process_args(int argc, char * const argv[], struct inputParams *p) {

    int next_option;
    const char* short_options = "c:hp:v";
	char paramfile[80];
    int modelSet=0, sourceSet=0, outputSet=0, no_option=1;
    do {
        next_option = getopt (argc, argv, short_options);
        switch (next_option)
        {
                
			case  'c' :
				strcpy(p->checkpointfile, optarg);
                p->checkpoint = 1;
                no_option = 0;
                break;
				
            case  'h' : /* -h or --help */
                /* User has requested usage information. Print it to standard
                 output, and exit with exit code zero (normal termination). */
                print_usage(argv[0], stdout, 0);
                
            case  'p':
                strcpy(paramfile, optarg);
                no_option = 0;
                break;
                
			case  'v':
				verbose++;
				break;
				
            case  '?' : /* The user specified an invalid option. */
                /* Print usage information to standard error, and exit with exit
                 code one (indicating abnormal termination). */
                print_usage (argv[0], stderr, 1);
                
            case -1: /* Done with options. */
                break;
                
            default: /* Something else: unexpected. */
                abort ();

        }
    } while (next_option != -1);
	
    if ( no_option==1 )
        print_usage (argv[0], stderr, 1);
	
	if ( p->checkpoint == 1 ) return;

	FILE *fid;
	fid = fopen(paramfile, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open parameter file %s for reading\n",
				paramfile);
        exit(1);
    }
	char line[200], value[50], keyword[50], *s1, *s2;
	while ( fgets( line, 200, fid) != NULL ) {
		if ( sscanf(line, "%s", value) <= 0 ) {
			continue;
		}		
		s1 = strchr(line, '#');
		s1++;
		s2 = strchr(line, ',');
		s2[0] = '\0';
		sscanf(s1, "%s", keyword);
		
		if ( strcmp(keyword, "model") == 0 ) {
			strcpy(p->modelfile, value);
			modelSet = 1;
		}
		else if ( strcmp(keyword, "source") == 0 ) {
			strcpy(p->sourcefile, value);
			sourceSet = 1;
		}
		else if ( strcmp(keyword, "output") == 0 ) {
			strcpy(p->outputfile, value);
			outputSet = 1;
		}
		else if ( strcmp(keyword, "basename") == 0 ) {
			strcpy(p->basename, value);
		}
		else if ( strcmp(keyword, "time") == 0 ) {
			p->duration = strtod(value, NULL);
		}
		else if ( strcmp(keyword, "dt") == 0 ) {
			p->dt = strtod(value, NULL) * 0.001;  // convert ms to s
		}
		else if ( strcmp(keyword, "abs") == 0 ) {
			p->iwipe = (short)atoi(value);
		}
		else if ( strcmp(keyword, "nl") == 0 ) {
			p->ab->np = (size_t)atoi(value);
		}
		else if ( strcmp(keyword, "plstr") == 0 ) {
			p->plotStrips = (short)atoi(value);
		}
		else if ( strcmp(keyword, "segy") == 0 ) {
			p->segy = (short)atoi(value);
		}
		else if ( strcmp(keyword, "check_model") == 0 ) {
			p->check_model = (short)atoi(value);
		}
		else if ( strcmp(keyword, "chkpt_inc") == 0 ) {
			p->chkpt_inc = atoi(value);
		}
		else if ( strcmp(keyword, "Rc") == 0 ) {
			p->ab->Rc = strtod(value, NULL);
		}
		else if ( strcmp(keyword, "kappa_max") == 0 ) {
			p->ab->kappa_max = strtod(value, NULL);
		}
		else if ( strcmp(keyword, "pml_median") == 0 ) {
			p->ab->median = (short)atoi(value);
		}
		else if ( strcmp(keyword, "alpha_max") == 0 ) {
			p->ab->alpha_max = (short)atoi(value);
		}
		else if ( strcmp(keyword, "saveEnergy") == 0 ) {
			p->saveEnergy = (short)atoi(value);
		}
		else if ( strcmp(keyword, "EnergyROI") == 0 ) {
			sscanf(value, "%lf/%lf/%lf/%lf", &(p->roiExmin), &(p->roiExmax),
				   &(p->roiEzmin), &(p->roiEzmax));
		}
		else if ( strcmp(keyword, "shotpoint") == 0 ) {
			p->shotpt_no = atoi(value);
		}
		else if ( strcmp(keyword, "simulateCMP") == 0 ) {
			p->simulateCMP = atoi(value);
		}
        else if ( strcmp(keyword, "azimuthal_mode") == 0 ) {
            p->n = atoi(value);
        }
        
	}
    if ( p->iwipe == 0 ) p->ab->np = 0;
}

void read_grid_params(const char filename[], struct grid *g) {
    FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    fscanf(fid, "%zd %zd", &(g->nx), &(g->nz));
    fscanf(fid, "%lf %lf", &(g->dx), &(g->dz));
    fscanf(fid, "%lf %lf", &(g->x0), &(g->z0));
    
	char type[30];
	fscanf(fid, "%s", type);
	for (size_t i=0; i<strlen(type); ++i ) type[i] = tolower(type[i]);
	//short with_pml=0;
	if ( strcmp(type, "isotropic_pml") == 0 ) {  // pml layers included in model file
		g->nx2 = g->nx;
		g->nx -= 2*g->ab.np;
		g->nz2 = g->nz;
		g->nz -= 2*g->ab.np;		
	}
    else if ( strcmp(type, "isotropic_elastic_cyl") == 0) {
        g->coord = CYLINDRICAL;
        g->nx2 = g->nx +   g->ab.np;
        g->nz2 = g->nz + 2*g->ab.np;
    }
	else {
		g->nx2 = g->nx + 2*g->ab.np;
		g->nz2 = g->nz + 2*g->ab.np;
	}
    fclose(fid);
}

void read_grid_params_ve_sh(const char filename[], struct grid *g,
							struct materialPropertiesVE_SH_VTI *mp) {
    FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    fscanf(fid, "%zd %zd", &(g->nx), &(g->nz));
    fscanf(fid, "%lf %lf", &(g->dx), &(g->dz));
    fscanf(fid, "%lf %lf", &(g->x0), &(g->z0));
    
	char type[30];
	fscanf(fid, "%s", type);
	for (size_t i=0; i<strlen(type); ++i ) type[i] = tolower(type[i]);
    if ( strcmp(type, "vti_sh_viscoelastic") != 0 ) {
        fprintf(stderr, "Error, model file %s not for a viscoelastic VTI medium\n", filename);
        exit(1);
    }
    
    g->nx2 = g->nx + 2*g->ab.np;
    g->nz2 = g->nz + 2*g->ab.np;
    fclose(fid);
}

void read_grid_params_ve(const char filename[], struct grid *g,
                         struct materialPropertiesVE_VTI *mp) {
    FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    fscanf(fid, "%zd %zd", &(g->nx), &(g->nz));
    fscanf(fid, "%lf %lf", &(g->dx), &(g->dz));
    fscanf(fid, "%lf %lf", &(g->x0), &(g->z0));
    
	char type[30];
	fscanf(fid, "%s", type);
	for (size_t i=0; i<strlen(type); ++i ) type[i] = tolower(type[i]);
    if ( strcmp(type, "vti_viscoelastic") == 0 ) {
        fscanf(fid, "%d", &(mp->L));
    } else {
        fprintf(stderr, "Error, model file %s not for a viscoelastic VTI medium\n", filename);
        exit(1);
    }
    
    g->nx2 = g->nx + 2*g->ab.np;
    g->nz2 = g->nz + 2*g->ab.np;
    fclose(fid);
}

void read_grid3d_params(const char filename[], struct grid3d *g) {
    FILE *fid;
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    
    fscanf(fid, "%zd %zd %zd", &(g->nx), &(g->ny), &(g->nz));
    fscanf(fid, "%lf %lf %lf", &(g->dx), &(g->dy), &(g->dz));
    fscanf(fid, "%lf %lf %lf", &(g->x0), &(g->y0), &(g->z0));
    
    g->nx2 = g->nx + 2*g->ab.np;
    g->ny2 = g->ny + 2*g->ab.np;
    g->nz2 = g->nz + 2*g->ab.np;
    fclose(fid);
    
}

void read_model(const char filename[], const struct grid *g,
                struct materialProperties *mp) {
    FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    fscanf(fid, "%*d %*d");
    fscanf(fid, "%*f %*f");
    fscanf(fid, "%*f %*f");
	char type[30];
	fscanf(fid, "%s", type);
	for (size_t i=0; i<strlen(type); ++i ) type[i] = tolower(type[i]);
	short with_pml=0;
	if ( strcmp(type, "isotropic_pml") == 0 ) {
		with_pml = 1;
	}
	else if ( strcmp(type, "isotropic") != 0 ) {
        fprintf(stderr, "Error, model in %s not isotropic\n", filename);
        exit(1);
    }
    
	int flag1D=0, nread;
	if ( with_pml == 0 ) {
		size_t nv=0;
		for (size_t i=0; i<g->nx; ++i) {
			size_t ind = (i+g->ab.np)*g->nz2 + g->ab.np;
			
			nread = fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
						   &(mp->K_m[ind]),
						   &(mp->K_s[ind]),
						   &(mp->K_f[ind]),
						   &(mp->phi[ind]),
						   &(mp->mu[ind]),
						   &(mp->rho_s[ind]),
						   &(mp->rho_f[ind]),
						   &(mp->T[ind]),
						   &(mp->eta[ind]),
						   &(mp->kappa[ind]),
						   &(mp->Q[ind]),
						   &(mp->f0[ind]));
			
			if ( nread<=0 && i==1 ) {
				flag1D = 1;
				break;
			}
			
			for (size_t j=1; j<g->nz; ++j, ++nv) {
				ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
				
				nread = fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
							   &(mp->K_m[ind]),
							   &(mp->K_s[ind]),
							   &(mp->K_f[ind]),
							   &(mp->phi[ind]),
							   &(mp->mu[ind]),
							   &(mp->rho_s[ind]),
							   &(mp->rho_f[ind]),
							   &(mp->T[ind]),
							   &(mp->eta[ind]),
							   &(mp->kappa[ind]),
							   &(mp->Q[ind]),
							   &(mp->f0[ind]));
				
				if (nread != 12) {
					fprintf(stderr, "Error: something wrong with file %s\n", filename);
					exit(1);
				}
			}
		}
		fclose(fid);
		
		if ( flag1D == 1 ) {
			for (size_t i=1; i<g->nx; ++i) {
				for (size_t j=0; j<g->nz; ++j, ++nv) {
					size_t ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
					size_t ind0 = g->ab.np*g->nz2 + j+g->ab.np;
					
					mp->K_m[ind]   = mp->K_m[ind0];
					mp->K_s[ind]   = mp->K_s[ind0];
					mp->K_f[ind]   = mp->K_f[ind0];
					mp->phi[ind]   = mp->phi[ind0];
					mp->mu[ind]    = mp->mu[ind0];
					mp->rho_s[ind] = mp->rho_s[ind0];
					mp->rho_f[ind] = mp->rho_f[ind0];
					mp->T[ind]     = mp->T[ind0];
					mp->eta[ind]   = mp->eta[ind0];
					mp->kappa[ind] = mp->kappa[ind0];
					mp->Q[ind]     = mp->Q[ind0];
					mp->f0[ind]    = mp->f0[ind0];
				}
			}
		}
		
		double Qend = 1.0;
		// Duplicate values at edges of the model inside the absorbing sides of the grid
		for (size_t i=0; i<g->nx; ++i) {
			for (size_t j=0; j<g->ab.np; ++j) {
				
				size_t i1 = (i+g->ab.np)*g->nz2 + j;
				size_t i2 = (i+g->ab.np)*g->nz2 + g->ab.np;
				
				mp->K_m[i1]   = mp->K_m[i2];
				mp->K_s[i1]   = mp->K_s[i2];
				mp->K_f[i1]   = mp->K_f[i2];
				mp->phi[i1]   = mp->phi[i2];
				mp->mu[i1]    = mp->mu[i2];
				mp->rho_s[i1] = mp->rho_s[i2];
				mp->rho_f[i1] = mp->rho_f[i2];
				mp->T[i1]     = mp->T[i2];
				mp->eta[i1]   = mp->eta[i2];
				mp->kappa[i1] = mp->kappa[i2];
				mp->f0[i1]    = mp->f0[i2];
				
				if ( g->Qdamping == 1 ) {
					double K = ( Qend - mp->Q[i1] ) / (g->ab.np*g->ab.np*g->ab.np);
					mp->Q[i1]     = mp->Q[i2] + K*(j+1)*(j+1)*(j+1);
				} else {
					mp->Q[i1]     = mp->Q[i2];
				}
				
				i1 = (i+g->ab.np)*g->nz2 + g->ab.np + g->nz + j;
				i2 = (i+g->ab.np)*g->nz2 + g->ab.np + g->nz - 1;
				
				mp->K_m[i1]   = mp->K_m[i2];
				mp->K_s[i1]   = mp->K_s[i2];
				mp->K_f[i1]   = mp->K_f[i2];
				mp->phi[i1]   = mp->phi[i2];
				mp->mu[i1]    = mp->mu[i2];
				mp->rho_s[i1] = mp->rho_s[i2];
				mp->rho_f[i1] = mp->rho_f[i2];
				mp->T[i1]     = mp->T[i2];
				mp->eta[i1]   = mp->eta[i2];
				mp->kappa[i1] = mp->kappa[i2];
				mp->f0[i1]    = mp->f0[i2];
				
				if ( g->Qdamping == 1 ) {
					double K = ( Qend - mp->Q[i1] ) / (g->ab.np*g->ab.np*g->ab.np);
					mp->Q[i1]     = mp->Q[i2] + K*(j+1)*(j+1)*(j+1);
				} else {
					mp->Q[i1]     = mp->Q[i2];
				}
			}
		}
		
		for (size_t i=0; i<g->ab.np; ++i) {
			for (size_t j=0; j<g->nz2; ++j) {
				
				size_t i1 = i*g->nz2 + j;
				size_t i2 = g->ab.np*g->nz2 + j;
				
				mp->K_m[i1]   = mp->K_m[i2];
				mp->K_s[i1]   = mp->K_s[i2];
				mp->K_f[i1]   = mp->K_f[i2];
				mp->phi[i1]   = mp->phi[i2];
				mp->mu[i1]    = mp->mu[i2];
				mp->rho_s[i1] = mp->rho_s[i2];
				mp->rho_f[i1] = mp->rho_f[i2];
				mp->T[i1]     = mp->T[i2];
				mp->eta[i1]   = mp->eta[i2];
				mp->kappa[i1] = mp->kappa[i2];
				mp->f0[i1]    = mp->f0[i2];
				
				if ( g->Qdamping == 1 ) {
					double K = ( Qend - mp->Q[i1] ) / (g->ab.np*g->ab.np*g->ab.np);
					mp->Q[i1]     = mp->Q[i2] + K*(j+1)*(j+1)*(j+1);
				} else {
					mp->Q[i1]     = mp->Q[i2];
				}
				
				i1 = (i+g->nx+g->ab.np)*g->nz2 + j;
				i2 = (g->nx+g->ab.np-1)*g->nz2 + j;
				
				mp->K_m[i1]   = mp->K_m[i2];
				mp->K_s[i1]   = mp->K_s[i2];
				mp->K_f[i1]   = mp->K_f[i2];
				mp->phi[i1]   = mp->phi[i2];
				mp->mu[i1]    = mp->mu[i2];
				mp->rho_s[i1] = mp->rho_s[i2];
				mp->rho_f[i1] = mp->rho_f[i2];
				mp->T[i1]     = mp->T[i2];
				mp->eta[i1]   = mp->eta[i2];
				mp->kappa[i1] = mp->kappa[i2];
				mp->f0[i1]    = mp->f0[i2];
				
				if ( g->Qdamping == 1 ) {
					double K = ( Qend - mp->Q[i1] ) / (g->ab.np*g->ab.np*g->ab.np);
					mp->Q[i1]     = mp->Q[i2] + K*(j+1)*(j+1)*(j+1);
				} else {
					mp->Q[i1]     = mp->Q[i2];
				}
			}
		}
	} else {
		for (size_t i=0, ind=0; i<g->nx2; ++i) {
			for (size_t j=0; j<g->nz2; ++j) {
				
				nread = fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
							   &(mp->K_m[ind]),
							   &(mp->K_s[ind]),
							   &(mp->K_f[ind]),
							   &(mp->phi[ind]),
							   &(mp->mu[ind]),
							   &(mp->rho_s[ind]),
							   &(mp->rho_f[ind]),
							   &(mp->T[ind]),
							   &(mp->eta[ind]),
							   &(mp->kappa[ind]),
							   &(mp->Q[ind]),
							   &(mp->f0[ind]));
				ind++;
				if (nread != 12 && ind<(g->nx2*g->nz2) ) {
					fprintf(stderr, "Error: something wrong with file, value %zd %s\n",
							ind, filename);
					exit(1);
				}
			}
		}
		fclose(fid);
	}
}

void read_modelVTI(const char filename[], const struct grid *g, 
                struct materialPropertiesVTI *mp) {
    FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    fscanf(fid, "%*d %*d");
    fscanf(fid, "%*f %*f");
    fscanf(fid, "%*f %*f");
	char type[30];
	fscanf(fid, "%s", type);
	for (size_t i=0; i<strlen(type); ++i ) type[i] = tolower(type[i]);
	short with_pml=0;
	if ( strcmp(type, "anisotropic__vti_pml") == 0 ) {
		with_pml = 1;
	}
	else if ( strcmp(type, "anisotropic_vti") != 0 ) {
        fprintf(stderr, "Error, model in %s not VTI anisotropic\n", filename);
        exit(1);
    }
    
    int flag1D=0, nread;
	double tmp1, tmp2;
	if ( with_pml == 0 ) {
		size_t nv=0;
		for (size_t i=0; i<g->nx; ++i) {
			size_t ind = (i+g->ab.np)*g->nz2 + g->ab.np;
        
			nread = fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
						   &(mp->K_m[ind]),
						   &(mp->K_s[ind]),
						   &(mp->K_f[ind]),
						   &(mp->phi[ind]),
						   &(mp->mu[ind]),
						   &(mp->rho_s[ind]),
						   &(mp->rho_f[ind]),
						   &(mp->T3[ind]),
						   &(mp->eta[ind]),
						   &(mp->kappa3[ind]),
						   &(mp->Q[ind]),
						   &(mp->f0[ind]),
						   &(mp->epsilon[ind]),
						   &(mp->delta[ind]),
						   &tmp1,
						   &tmp2);
			
			mp->kappa1[ind] = mp->kappa3[ind] / tmp1;
			mp->T1[ind] = mp->T3[ind] / tmp2;
			
			if ( nread<=0 && i==1 ) {
				flag1D = 1;
				break;
			}
			
			for (size_t j=1; j<g->nz; ++j, ++nv) {
				ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
				
				nread = fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
							   &(mp->K_m[ind]),
							   &(mp->K_s[ind]),
							   &(mp->K_f[ind]),
							   &(mp->phi[ind]),
							   &(mp->mu[ind]),
							   &(mp->rho_s[ind]),
							   &(mp->rho_f[ind]),
							   &(mp->T3[ind]),
							   &(mp->eta[ind]),
							   &(mp->kappa3[ind]),
							   &(mp->Q[ind]),
							   &(mp->f0[ind]),
							   &(mp->epsilon[ind]),
							   &(mp->delta[ind]),
							   &tmp1,
							   &tmp2);
				
				mp->kappa1[ind] = mp->kappa3[ind] / tmp1;
				mp->T1[ind] = mp->T3[ind] / tmp2;
				
				if ( nread != 16 ) {
					fprintf(stderr, "Error: something wrong with file %s\n", filename);
					exit(1);
				}
			}
		}
		fclose(fid);
		
		if ( flag1D == 1 ) {
			for (size_t i=1; i<g->nx; ++i) {
				for (size_t j=0; j<g->nz; ++j, ++nv) {
					size_t ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
					size_t ind0 = g->ab.np*g->nz2 + j+g->ab.np;
					
					mp->K_m[ind]     = mp->K_m[ind0];
					mp->K_s[ind]     = mp->K_s[ind0];
					mp->K_f[ind]     = mp->K_f[ind0];
					mp->phi[ind]     = mp->phi[ind0];
					mp->mu[ind]      = mp->mu[ind0];
					mp->rho_s[ind]   = mp->rho_s[ind0];
					mp->rho_f[ind]   = mp->rho_f[ind0];
					mp->T1[ind]      = mp->T1[ind0];
					mp->T3[ind]      = mp->T3[ind0];
					mp->eta[ind]     = mp->eta[ind0];
					mp->kappa1[ind]  = mp->kappa1[ind0];
					mp->kappa3[ind]  = mp->kappa3[ind0];
					mp->Q[ind]       = mp->Q[ind0];
					mp->epsilon[ind] = mp->epsilon[ind0];
					mp->delta[ind]   = mp->delta[ind0];
					mp->f0[ind]      = mp->f0[ind0];                
				}
			}
		}
		
		double Qend = 1.0;
		// Duplicate values at edges of the model inside the absorbing sides of the grid
		for (size_t i=0; i<g->nx; ++i) {
			for (size_t j=0; j<g->ab.np; ++j) {
				
				size_t i1 = (i+g->ab.np)*g->nz2 + j;
				size_t i2 = (i+g->ab.np)*g->nz2 + g->ab.np;
				
				mp->K_m[i1]     = mp->K_m[i2];
				mp->K_s[i1]     = mp->K_s[i2];
				mp->K_f[i1]     = mp->K_f[i2];
				mp->phi[i1]     = mp->phi[i2];
				mp->mu[i1]      = mp->mu[i2];
				mp->rho_s[i1]   = mp->rho_s[i2];
				mp->rho_f[i1]   = mp->rho_f[i2];
				mp->T1[i1]      = mp->T1[i2];
				mp->T3[i1]      = mp->T3[i2];
				mp->eta[i1]     = mp->eta[i2];
				mp->kappa1[i1]  = mp->kappa1[i2];
				mp->kappa3[i1]  = mp->kappa3[i2];
				mp->epsilon[i1] = mp->epsilon[i2];
				mp->delta[i1]   = mp->delta[i2];
				mp->f0[i1]      = mp->f0[i2];
				
				if ( g->Qdamping == 1 ) {
					double K = ( Qend - mp->Q[i1] ) / (g->ab.np*g->ab.np*g->ab.np);
					mp->Q[i1] = mp->Q[i2] + K*(j+1)*(j+1)*(j+1);
				} else {
					mp->Q[i1] = mp->Q[i2];
				}
				
				i1 = (i+g->ab.np)*g->nz2 + g->ab.np + g->nz + j;
				i2 = (i+g->ab.np)*g->nz2 + g->ab.np + g->nz - 1;
				
				mp->K_m[i1]     = mp->K_m[i2];
				mp->K_s[i1]     = mp->K_s[i2];
				mp->K_f[i1]     = mp->K_f[i2];
				mp->phi[i1]     = mp->phi[i2];
				mp->mu[i1]      = mp->mu[i2];
				mp->rho_s[i1]   = mp->rho_s[i2];
				mp->rho_f[i1]   = mp->rho_f[i2];
				mp->T1[i1]      = mp->T1[i2];
				mp->T3[i1]      = mp->T3[i2];
				mp->eta[i1]     = mp->eta[i2];
				mp->kappa1[i1]  = mp->kappa1[i2];
				mp->kappa3[i1]  = mp->kappa3[i2];
				mp->epsilon[i1] = mp->epsilon[i2];
				mp->delta[i1]   = mp->delta[i2];
				mp->f0[i1]      = mp->f0[i2];
				
				if ( g->Qdamping == 1 ) {
					double K = ( Qend - mp->Q[i1] ) / (g->ab.np*g->ab.np*g->ab.np);
					mp->Q[i1] = mp->Q[i2] + K*(j+1)*(j+1)*(j+1);
				} else {
					mp->Q[i1] = mp->Q[i2];
				}
			}
		}
    
		for (size_t i=0; i<g->ab.np; ++i) {
			for (size_t j=0; j<g->nz2; ++j) {
				
				size_t i1 = i*g->nz2 + j;
				size_t i2 = g->ab.np*g->nz2 + j;
				
				mp->K_m[i1]     = mp->K_m[i2];
				mp->K_s[i1]     = mp->K_s[i2];
				mp->K_f[i1]     = mp->K_f[i2];
				mp->phi[i1]     = mp->phi[i2];
				mp->mu[i1]      = mp->mu[i2];
				mp->rho_s[i1]   = mp->rho_s[i2];
				mp->rho_f[i1]   = mp->rho_f[i2];
				mp->T1[i1]      = mp->T1[i2];
				mp->T3[i1]      = mp->T3[i2];
				mp->eta[i1]     = mp->eta[i2];
				mp->kappa1[i1]  = mp->kappa1[i2];
				mp->kappa3[i1]  = mp->kappa3[i2];
				mp->epsilon[i1] = mp->epsilon[i2];
				mp->delta[i1]   = mp->delta[i2];
				mp->f0[i1]      = mp->f0[i2];
				
				if ( g->Qdamping == 1 ) {
					double K = ( Qend - mp->Q[i1] ) / (g->ab.np*g->ab.np*g->ab.np);
					mp->Q[i1] = mp->Q[i2] + K*(j+1)*(j+1)*(j+1);
				} else {
					mp->Q[i1] = mp->Q[i2];
				}
				
				i1 = (i+g->nx+g->ab.np)*g->nz2 + j;
				i2 = (g->nx+g->ab.np-1)*g->nz2 + j;
				
				mp->K_m[i1]     = mp->K_m[i2];
				mp->K_s[i1]     = mp->K_s[i2];
				mp->K_f[i1]     = mp->K_f[i2];
				mp->phi[i1]     = mp->phi[i2];
				mp->mu[i1]      = mp->mu[i2];
				mp->rho_s[i1]   = mp->rho_s[i2];
				mp->rho_f[i1]   = mp->rho_f[i2];
				mp->T1[i1]      = mp->T1[i2];
				mp->T3[i1]      = mp->T3[i2];
				mp->eta[i1]     = mp->eta[i2];
				mp->kappa1[i1]  = mp->kappa1[i2];
				mp->kappa3[i1]  = mp->kappa3[i2];
				mp->epsilon[i1] = mp->epsilon[i2];
				mp->delta[i1]   = mp->delta[i2];
				mp->f0[i1]      = mp->f0[i2];
				
				if ( g->Qdamping == 1 ) {
					double K = ( Qend - mp->Q[i1] ) / (g->ab.np*g->ab.np*g->ab.np);
					mp->Q[i1] = mp->Q[i2] + K*(j+1)*(j+1)*(j+1);
				} else {
					mp->Q[i1] = mp->Q[i2];
				}
			}
		}
	} else {
		for (size_t i=0, ind=0; i<g->nx2; ++i) {
			for (size_t j=0; j<g->nz2; ++j, ind++) {
				
				nread = fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
							   &(mp->K_m[ind]),
							   &(mp->K_s[ind]),
							   &(mp->K_f[ind]),
							   &(mp->phi[ind]),
							   &(mp->mu[ind]),
							   &(mp->rho_s[ind]),
							   &(mp->rho_f[ind]),
							   &(mp->T3[ind]),
							   &(mp->eta[ind]),
							   &(mp->kappa3[ind]),
							   &(mp->Q[ind]),
							   &(mp->f0[ind]),
							   &(mp->epsilon[ind]),
							   &(mp->delta[ind]),
							   &tmp1,
							   &tmp2);
				
				mp->kappa1[ind] = mp->kappa3[ind] / tmp1;
				mp->T1[ind] = mp->T3[ind] / tmp2;
				
				if ( nread != 16 ) {
					fprintf(stderr, "Error: something wrong with file %s\n", filename);
					exit(1);
				}				
			}
		}
		fclose(fid);
	}
}
	
void read_model_e(const char filename[], const struct grid *g, double *vp,
				  double *vs, double *rho) {
    FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    fscanf(fid, "%*d %*d");
    fscanf(fid, "%*f %*f");
    fscanf(fid, "%*f %*f");
	char type[30];
	fscanf(fid, "%s", type);
	if ( strcmp(type, "isotropic_elastic") != 0 ) {
        fprintf(stderr, "Error, model in %s not isotropic_elastic\n", filename);
        exit(1);
    }
    
    size_t nv=0;
    int flag1D=0, nread;
    for (size_t i=0; i<g->nx; ++i) {
        size_t ind = i*g->nz;
        
        nread = fscanf(fid, "%lf %lf %lf", &(vp[ind]), &(vs[ind]), &(rho[ind]));
        
        if ( nread<=0 && i==1 ) {
            flag1D = 1;
            break;
        }
        
        for (size_t j=1; j<g->nz; ++j, ++nv) {
            ind = i*g->nz + j;
			
            nread = fscanf(fid, "%lf %lf %lf", &(vp[ind]), &(vs[ind]), &(rho[ind]));
			
			if (nread != 3) {
				fprintf(stderr, "Error: something wrong with file %s\n", filename);
				exit(1);
			}
        }
    }
    fclose(fid);
    
    if ( flag1D == 1 ) {
        for (size_t i=1; i<g->nx; ++i) {
            for (size_t j=0; j<g->nz; ++j, ++nv) {
                size_t ind = i*g->nz + j;
				
                vp[ind]   = vp[j];
                vs[ind]   = vs[j];
                rho[ind]  = rho[j];				
            }
        }
    }
}

void read_model_e_cyl(const char filename[], const struct grid *g, double *vp,
                      double *vs, double *rho) {
    FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    fscanf(fid, "%*d %*d");
    fscanf(fid, "%*f %*f");
    fscanf(fid, "%*f %*f");
    char type[30];
    fscanf(fid, "%s", type);
    if ( strcmp(type, "isotropic_elastic_cyl") != 0 ) {
        fprintf(stderr, "Error, model in %s not isotropic_elastic_cyl\n", filename);
        exit(1);
    }
    
    size_t nv=0;
    int flag1D=0, nread;
    size_t i=0;
    for (size_t j=0; j<g->nz; ++j) {
        size_t ind = i*g->nz2 + j+g->ab.np;
        
        nread = fscanf(fid, "%lf %lf %lf", &(vp[ind]), &(vs[ind]), &(rho[ind]));
        
        if ( nread<=0 && j==1 ) {
            flag1D = 1;
            break;
        }
        
        for (i=1; i<g->nx; ++i, ++nv) {
            ind = i*g->nz2 + j+g->ab.np;
            
            nread = fscanf(fid, "%lf %lf %lf", &(vp[ind]), &(vs[ind]), &(rho[ind]));
            
            if (nread != 3) {
                fprintf(stderr, "Error: something wrong with file %s\n", filename);
                exit(1);
            }
        }
        i = 0;
    }
    fclose(fid);
    
    if ( flag1D == 1 ) {
        for (size_t i=0; i<g->nx; ++i) {
            for (size_t j=1; j<g->nz; ++j, ++nv) {
                size_t ind = i*g->nz2 + j+g->ab.np;
                
                vp[ind]   = vp[i*g->nz2  + g->ab.np];
                vs[ind]   = vs[i*g->nz2  + g->ab.np];
                rho[ind]  = rho[i*g->nz2 + g->ab.np];
            }
        }
    }
    
    // duplicate values in PML
    for ( size_t i=0; i<g->nx; ++i ) {
        for ( size_t j=0; j<g->ab.np; ++j ) {
            vp[i*g->nz2+j]  = vp[i*g->nz2  + g->ab.np];
            vs[i*g->nz2+j]  = vs[i*g->nz2  + g->ab.np];
            rho[i*g->nz2+j] = rho[i*g->nz2 + g->ab.np];
        }
        for ( size_t j=g->nz+g->ab.np; j<g->nz2; ++j ) {
            vp[i*g->nz2+j]  = vp[i*g->nz2  + g->nz+g->ab.np-1];
            vs[i*g->nz2+j]  = vs[i*g->nz2  + g->nz+g->ab.np-1];
            rho[i*g->nz2+j] = rho[i*g->nz2 + g->nz+g->ab.np-1];
        }
    }
    for ( size_t i=g->nx; i<g->nx2; ++i ) {
        for ( size_t j=0; j<g->nz2; ++j ) {
            vp[i*g->nz2+j]  = vp[(g->nx-1)*g->nz2  + j];
            vs[i*g->nz2+j]  = vs[(g->nx-1)*g->nz2  + j];
            rho[i*g->nz2+j] = rho[(g->nx-1)*g->nz2 + j];
        }
    }
}

void read_model_a(const char filename[], const struct grid3d *g, double *vp,
				  double *rho) {
    FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    fscanf(fid, "%*d %*d %*d");
    fscanf(fid, "%*f %*f %*f");
    fscanf(fid, "%*f %*f %*f");
	char type[30];
	fscanf(fid, "%s", type);
	if ( strcmp(type, "isotropic_acoustic_3d") != 0 ) {
        fprintf(stderr, "Error, model in %s not isotropic_acoustic_3d\n", filename);
        exit(1);
    }
    
    int flag1D=0, nread;
    for (size_t i=0; i<g->nx; ++i) {
        for (size_t j=0; j<g->ny; ++j) {

            size_t ind = ((i+g->ab.np)*g->ny2+j+g->ab.np)*g->nz2+g->ab.np;
        
            nread = fscanf(fid, "%lf %lf", &(vp[ind]), &(rho[ind]));
        
            if ( nread<=0 && j==1 ) {
                flag1D = 1;
                break;
            }
        
            for (size_t k=1; k<g->nz; ++k) {
                ind = ((i+g->ab.np)*g->ny2+j+g->ab.np)*g->nz2 + k+g->ab.np;
                
                nread = fscanf(fid, "%lf %lf", &(vp[ind]), &(rho[ind]));
                
                if (nread != 2) {
                    fprintf(stderr, "Error: something wrong with file %s\n", filename);
                    exit(1);
                }
            }
        }
        if ( flag1D ) break;
    }
    fclose(fid);
    
    if ( flag1D == 1 ) {
        for (size_t i=0; i<g->nx; ++i) {
            for (size_t j=0; j<g->ny; ++j) {
                for (size_t k=0; k<g->nz; ++k) {
                    size_t ind = ((i+g->ab.np)*g->ny2+j+g->ab.np)*g->nz2 + k+g->ab.np;
				
                    vp[ind]   =  vp[(g->ab.np*g->ny2+g->ab.np)*g->nz2 + k+g->ab.np];
                    rho[ind]  = rho[(g->ab.np*g->ny2+g->ab.np)*g->nz2 + k+g->ab.np];
                }
            }
        }
    }
    
    // Duplicate values at edges of the model inside the absorbing sides of the grid
    for (size_t i=0; i<g->nx; ++i) {
        for (size_t j=0; j<g->ny; ++j) {
            for (size_t k=0; k<g->ab.np; ++k) {
            
                size_t i1 = ((i+g->ab.np)*g->ny2+j+g->ab.np)*g->nz2 + k;
                size_t i2 = ((i+g->ab.np)*g->ny2+j+g->ab.np)*g->nz2 + g->ab.np;
                
                vp[i1] =  vp[i2];
                rho[i1] = rho[i2];
                
                i1 = ((i+g->ab.np)*g->ny2+j+g->ab.np)*g->nz2 + g->ab.np + g->nz + k;
                i2 = ((i+g->ab.np)*g->ny2+j+g->ab.np)*g->nz2 + g->ab.np + g->nz - 1;
                
                vp[i1] =  vp[i2];
                rho[i1] = rho[i2];
            }
        }
    }
    for (size_t i=0; i<g->nx; ++i) {
        for (size_t k=0; k<g->nz2; ++k) {
            for ( size_t j=0; j<g->ab.np; ++j ) {
                
                size_t i1 = ((i+g->ab.np)*g->ny2+j)*g->nz2 + k;
                size_t i2 = ((i+g->ab.np)*g->ny2+g->ab.np)*g->nz2 + k;

                vp[i1] =  vp[i2];
                rho[i1] = rho[i2];
                
                i1 = ((i+g->ab.np)*g->ny2+g->ab.np+g->ny+j)*g->nz2 + k;
                i2 = ((i+g->ab.np)*g->ny2+g->ab.np+g->ny-1)*g->nz2 + k;

                vp[i1] =  vp[i2];
                rho[i1] = rho[i2];
            }
        }
    }
    for (size_t j=0; j<g->ny2; ++j) {
        for (size_t k=0; k<g->nz2; ++k) {
            for ( size_t i=0; i<g->ab.np; ++i ) {
                
                size_t i1 = (i*g->ny2+j)*g->nz2 + k;
                size_t i2 = (g->ab.np*g->ny2+j)*g->nz2 + k;

                vp[i1] =  vp[i2];
                rho[i1] = rho[i2];
                
                i1 = ((g->ab.np+g->nx+i)*g->ny2+j)*g->nz2 + k;
                i2 = ((g->ab.np+g->nx-1)*g->ny2+j)*g->nz2 + k;
                
                vp[i1] =  vp[i2];
                rho[i1] = rho[i2];
            }
        }
    }
	
//	for ( size_t i=0, n=0; i<g->nx2; ++i ) {
//		for ( size_t j=0; j<g->ny2; ++j ) {
//			for ( size_t k=0; k<g->nz2; ++k, ++n ) {
//				printf("%zd  %zd  %zd  %lf  %lf\n", i, j, k, vp[n], rho[n]);
//			}
//		}
//	}
	
}

void read_model_ve(const char filename[], const struct grid *g,
				   struct materialPropertiesVE_VTI *mp) {
    FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    fscanf(fid, "%*d %*d");
    fscanf(fid, "%*f %*f");
    fscanf(fid, "%*f %*f");
	char type[30];
	fscanf(fid, "%s", type);
	if ( strcmp(type, "vti_viscoelastic") != 0 ) {
        fprintf(stderr, "Error, model in %s not vti_viscoelastic\n", filename);
        exit(1);
    }
    fscanf(fid, "%*d");
    
    int flag1D=0, nread;
    for (size_t i=0; i<g->nx; ++i) {
        
        size_t ind = (i+g->ab.np)*g->nz2 + g->ab.np;
        
        double Vp0;
        double Vs0;
        
        nread = fscanf(fid, "%lf %lf %lf %lf %lf", &Vp0,
               &Vs0, &(mp->rho[ind]), &(mp->epsilon[ind]),
               &(mp->delta[ind]));
        
        mp->mu[ind] = Vs0*Vs0*mp->rho[ind];  // in Pa
        mp->K[ind] = 1.e-9 * (Vp0*Vp0*mp->rho[ind] - 4./3. * mp->mu[ind]);
        mp->mu[ind] *= 1.e-9;  // now in GPa
        
        if ( nread<=0 && i==1 ) {
            flag1D = 1;
            break;
        }

        for ( size_t nm=0; nm<mp->L; ++nm ) {
            fscanf(fid, "%lf %lf", &(mp->Q1[nm][ind]), &(mp->f1[nm][ind]));
        }
        fscanf(fid, "%lf %lf", &(mp->Q2[ind]), &(mp->f2[ind]));
        
        for (size_t j=1; j<g->nz; ++j) {
            ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
            fscanf(fid, "%lf %lf %lf %lf %lf", &Vp0,
                   &Vs0, &(mp->rho[ind]), &(mp->epsilon[ind]),
                   &(mp->delta[ind]));
            
            mp->mu[ind] = Vs0*Vs0*mp->rho[ind];  // in Pa
            mp->K[ind] = 1.e-9 * (Vp0*Vp0*mp->rho[ind] - 4./3. * mp->mu[ind]);
            mp->mu[ind] *= 1.e-9;  // now in GPa
            
            for ( size_t nm=0; nm<mp->L; ++nm ) {
                fscanf(fid, "%lf %lf", &(mp->Q1[nm][ind]), &(mp->f1[nm][ind]));
            }
            fscanf(fid, "%lf %lf", &(mp->Q2[ind]), &(mp->f2[ind]));
        }

        if ( flag1D ) break;
    }
    fclose(fid);
    
    if ( flag1D == 1 ) {
        for (size_t i=1; i<g->nx; ++i) {
            for (size_t j=0; j<g->nz; ++j) {
                size_t ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
                size_t ind0 = g->ab.np*g->nz2 + j+g->ab.np;
                
                mp->K[ind]       = mp->K[ind0];
                mp->mu[ind]      = mp->mu[ind0];
                mp->rho[ind]     = mp->rho[ind0];
                mp->epsilon[ind] = mp->epsilon[ind0];
                mp->delta[ind]   = mp->delta[ind0];
                for ( size_t nm=0; nm<mp->L; ++nm ) {
                    mp->Q1[nm][ind]      = mp->Q1[nm][ind0];
                    mp->f1[nm][ind]      = mp->f1[nm][ind0];
                }
                mp->Q2[ind]      = mp->Q2[ind0];
                mp->f2[ind]      = mp->f2[ind0];
            }
        }
    }
    
    // Duplicate values at edges of the model inside the absorbing sides of the grid
    for (size_t i=0; i<g->nx; ++i) {
        for (size_t j=0; j<g->ab.np; ++j) {
            
            size_t i1 = (i+g->ab.np)*g->nz2 + j;
            size_t i2 = (i+g->ab.np)*g->nz2 + g->ab.np;
            
            mp->K[i1]       = mp->K[i2];
            mp->mu[i1]      = mp->mu[i2];
            mp->rho[i1]     = mp->rho[i2];
            mp->epsilon[i1] = mp->epsilon[i2];
            mp->delta[i1]   = mp->delta[i2];
            for ( size_t nm=0; nm<mp->L; ++nm ) {
                mp->Q1[nm][i1]      = mp->Q1[nm][i2];
                mp->f1[nm][i1]      = mp->f1[nm][i2];
            }
            mp->Q2[i1]      = mp->Q2[i2];
            mp->f2[i1]      = mp->f2[i2];
            
            i1 = (i+g->ab.np)*g->nz2 + g->ab.np + g->nz + j;
            i2 = (i+g->ab.np)*g->nz2 + g->ab.np + g->nz - 1;
            
            mp->K[i1]       = mp->K[i2];
            mp->mu[i1]      = mp->mu[i2];
            mp->rho[i1]     = mp->rho[i2];
            mp->epsilon[i1] = mp->epsilon[i2];
            mp->delta[i1]   = mp->delta[i2];
            for ( size_t nm=0; nm<mp->L; ++nm ) {
                mp->Q1[nm][i1]      = mp->Q1[nm][i2];
                mp->f1[nm][i1]      = mp->f1[nm][i2];
            }
            mp->Q2[i1]      = mp->Q2[i2];
            mp->f2[i1]      = mp->f2[i2];
            
        }
    }
    
    for (size_t i=0; i<g->ab.np; ++i) {
        for (size_t j=0; j<g->nz2; ++j) {
            
            size_t i1 = i*g->nz2 + j;
            size_t i2 = g->ab.np*g->nz2 + j;
            
            mp->K[i1]       = mp->K[i2];
            mp->mu[i1]      = mp->mu[i2];
            mp->rho[i1]     = mp->rho[i2];
            mp->epsilon[i1] = mp->epsilon[i2];
            mp->delta[i1]   = mp->delta[i2];
            for ( size_t nm=0; nm<mp->L; ++nm ) {
                mp->Q1[nm][i1]      = mp->Q1[nm][i2];
                mp->f1[nm][i1]      = mp->f1[nm][i2];
            }
            mp->Q2[i1]      = mp->Q2[i2];
            mp->f2[i1]      = mp->f2[i2];
            
            i1 = (i+g->nx+g->ab.np)*g->nz2 + j;
            i2 = (g->nx+g->ab.np-1)*g->nz2 + j;
            
            mp->K[i1]       = mp->K[i2];
            mp->mu[i1]      = mp->mu[i2];
            mp->rho[i1]     = mp->rho[i2];
            mp->epsilon[i1] = mp->epsilon[i2];
            mp->delta[i1]   = mp->delta[i2];
            for ( size_t nm=0; nm<mp->L; ++nm ) {
                mp->Q1[nm][i1]      = mp->Q1[nm][i2];
                mp->f1[nm][i1]      = mp->f1[nm][i2];
            }
            mp->Q2[i1]      = mp->Q2[i2];
            mp->f2[i1]      = mp->f2[i2];
            
        }
    }
}

void read_model_ve_sh(const char filename[], const struct grid *g,
					  struct materialPropertiesVE_SH_VTI *mp) {
    FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    fscanf(fid, "%*d %*d");
    fscanf(fid, "%*f %*f");
    fscanf(fid, "%*f %*f");
	char type[30];
	fscanf(fid, "%s", type);
	if ( strcmp(type, "vti_sh_viscoelastic") != 0 ) {
        fprintf(stderr, "Error, model in %s not vti_sh_viscoelastic\n", filename);
        exit(1);
    }
    
    int flag1D=0, nread;
    for (size_t i=0; i<g->nx; ++i) {
        
        size_t ind = (i+g->ab.np)*g->nz2 + g->ab.np;
        
        double Vs0;
        
        nread = fscanf(fid, "%lf %lf %lf %lf %lf", &Vs0,
					   &(mp->rho[ind]), &(mp->gamma[ind]),
					   &(mp->Q[ind]), &(mp->f[ind]));

        mp->mu[ind] = Vs0*Vs0*mp->rho[ind] * 1.e-9;
        
		if ( nread<=0 && i==1 ) {
            flag1D = 1;
            break;
        }

        for (size_t j=1; j<g->nz; ++j) {
            ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
            fscanf(fid, "%lf %lf %lf %lf %lf", 
                   &Vs0, &(mp->rho[ind]), &(mp->gamma[ind]),
                   &(mp->Q[ind]), &(mp->f[ind]));
        
            mp->mu[ind] = Vs0*Vs0*mp->rho[ind] * 1.e-9;
        }
		
        if ( flag1D ) break;
    }
    fclose(fid);
    
    if ( flag1D == 1 ) {
        for (size_t i=1; i<g->nx; ++i) {
            for (size_t j=0; j<g->nz; ++j) {
                size_t ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
                size_t ind0 = g->ab.np*g->nz2 + j+g->ab.np;
                
                mp->mu[ind]      = mp->mu[ind0];
                mp->rho[ind]     = mp->rho[ind0];
                mp->gamma[ind]   = mp->gamma[ind0];
                mp->Q[ind]       = mp->Q[ind0];
                mp->f[ind]       = mp->f[ind0];
            }
        }
    }
    
    // Duplicate values at edges of the model inside the absorbing sides of the grid
    for (size_t i=0; i<g->nx; ++i) {
        for (size_t j=0; j<g->ab.np; ++j) {
            
            size_t i1 = (i+g->ab.np)*g->nz2 + j;
            size_t i2 = (i+g->ab.np)*g->nz2 + g->ab.np;
            
            mp->mu[i1]      = mp->mu[i2];
            mp->rho[i1]     = mp->rho[i2];
            mp->gamma[i1]   = mp->gamma[i2];
            mp->Q[i1]       = mp->Q[i2];
            mp->f[i1]       = mp->f[i2];
            
            i1 = (i+g->ab.np)*g->nz2 + g->ab.np + g->nz + j;
            i2 = (i+g->ab.np)*g->nz2 + g->ab.np + g->nz - 1;
            
            mp->mu[i1]      = mp->mu[i2];
            mp->rho[i1]     = mp->rho[i2];
            mp->gamma[i1]   = mp->gamma[i2];
            mp->Q[i1]       = mp->Q[i2];
            mp->f[i1]       = mp->f[i2];
            
        }
    }
    
    for (size_t i=0; i<g->ab.np; ++i) {
        for (size_t j=0; j<g->nz2; ++j) {
            
            size_t i1 = i*g->nz2 + j;
            size_t i2 = g->ab.np*g->nz2 + j;
            
            mp->mu[i1]      = mp->mu[i2];
            mp->rho[i1]     = mp->rho[i2];
            mp->gamma[i1]   = mp->gamma[i2];
            mp->Q[i1]       = mp->Q[i2];
            mp->f[i1]       = mp->f[i2];
            
            i1 = (i+g->nx+g->ab.np)*g->nz2 + j;
            i2 = (g->nx+g->ab.np-1)*g->nz2 + j;
            
            mp->mu[i1]      = mp->mu[i2];
            mp->rho[i1]     = mp->rho[i2];
            mp->gamma[i1]   = mp->gamma[i2];
            mp->Q[i1]       = mp->Q[i2];
            mp->f[i1]       = mp->f[i2];
            
        }
    }
}


void read_source(const char filename[], const struct grid *g, struct sourceParams *src) {
	FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }

	fscanf(fid, "%zd", &(src->nsrc));
	int nread = fscanf(fid, "%zd", &(src->nTemplate));
    if ( nread == 0 ) {
        src->nTemplate = 1;
    }
	if ( NULL == ( src->s = (struct source *) malloc(src->nsrc*src->nTemplate*sizeof(struct source)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }

	char t[50];
	for ( size_t n=0; n<src->nsrc*src->nTemplate; n+=src->nTemplate ) {
		
		fscanf(fid, "%s", t);
		if ( strcmp(t,"sx")==0 || strcmp(t,"SX")==0 || strcmp(t,"Sx")==0 ) {
			src->s[n].type = SX;
		} else if ( strcmp(t,"sz")==0 || strcmp(t,"SZ")==0 || strcmp(t,"Sz")==0 ) {
			src->s[n].type = SZ;
		} else if ( strcmp(t,"sxy")==0 || strcmp(t,"SXY")==0 || strcmp(t,"Sxy")==0 ) {
			src->s[n].type = SXY;
		} else if ( strcmp(t,"sxz")==0 || strcmp(t,"SXZ")==0 || strcmp(t,"Sxz")==0 ) {
			src->s[n].type = SXZ;
        } else if ( strcmp(t,"srz")==0 || strcmp(t,"SRZ")==0 || strcmp(t,"Srz")==0 ) {
            src->s[n].type = SXZ;
		} else if ( strcmp(t,"syz")==0 || strcmp(t,"SYZ")==0 || strcmp(t,"Syz")==0 ) {
			src->s[n].type = SYZ;
		} else if ( strcmp(t,"sf")==0 || strcmp(t,"SF")==0 || strcmp(t,"Sf")==0 ) {
			src->s[n].type = SF;
		} else if ( strcmp(t,"bulk")==0 || strcmp(t,"BULK")==0 || strcmp(t,"Bulk")==0 ) {
			src->s[n].type = BULK;
		} else if ( strcmp(t,"bulk_s")==0 || strcmp(t,"BULK_S")==0 || strcmp(t,"Bulk_s")==0 ) {
			src->s[n].type = BULK_S;
        } else if( strcmp(t,"fr")==0 || strcmp(t,"FR")==0 || strcmp(t,"Fr")==0 ) {
            src->s[n].type = FR;
        } else if( strcmp(t,"fz")==0 || strcmp(t,"FZ")==0 || strcmp(t,"Fz")==0 ) {
            src->s[n].type = FZ;
        } else if( strcmp(t,"ftheta")==0 || strcmp(t,"FTHETA")==0 || strcmp(t,"Ftheta")==0 ) {
            src->s[n].type = FT;
        } else if ( strcmp(t,"kurkjian")==0 || strcmp(t,"Kurkjian")==0 ) {
            src->s[n].type = KURKJIAN;
        } else {
            fprintf(stderr, "Error: type of source (%s) not implemented\n", t);
            exit(1);
        }
		
		fscanf(fid,"%lf", &(src->s[n].A) );
		fscanf(fid,"%lf", &(src->s[n].f) );
		fscanf(fid,"%lf %lf", &(src->s[n].x), &(src->s[n].z) );
		
//		src->s[n].A *= 1.e6;   // source strength now in Pa
        size_t i = round((src->s[n].x-g->x0)/g->dx)+g->ab.np;
        if ( g->coord == CYLINDRICAL )
            i = round((src->s[n].x-g->x0)/g->dx);
		size_t j = round((src->s[n].z-g->z0)/g->dz)+g->ab.np;
		src->s[n].i = i*g->nz2+j;
		
		if ( src->nTemplate == 9 ) {
			double A = src->s[n].A;
			src->s[n].A *= 0.5;
			src->s[n+1].A = 0.125*A;
			src->s[n+2].A = 0.25*A;
			src->s[n+3].A = 0.125*A;
			src->s[n+4].A = 0.25*A;
			src->s[n+5].A = 0.25*A;
			src->s[n+6].A = 0.125*A;
			src->s[n+7].A = 0.25*A;
			src->s[n+8].A = 0.125*A;

            src->s[n+1].f = src->s[n].f;
            src->s[n+2].f = src->s[n].f;
            src->s[n+3].f = src->s[n].f;
            src->s[n+4].f = src->s[n].f;
            src->s[n+5].f = src->s[n].f;
            src->s[n+6].f = src->s[n].f;
            src->s[n+7].f = src->s[n].f;
            src->s[n+8].f = src->s[n].f;

			i = round((src->s[n].x-g->x0-g->dx)/g->dx)+g->ab.np;
			j = round((src->s[n].z-g->z0-g->dz)/g->dz)+g->ab.np;
			src->s[n+1].i = i*g->nz2+j;
			i = round((src->s[n].x-g->x0      )/g->dx)+g->ab.np;
			src->s[n+2].i = i*g->nz2+j;
			i = round((src->s[n].x-g->x0+g->dx)/g->dx)+g->ab.np;
			src->s[n+3].i = i*g->nz2+j;
			i = round((src->s[n].x-g->x0-g->dx)/g->dx)+g->ab.np;
			j = round((src->s[n].z-g->z0      )/g->dz)+g->ab.np;
			src->s[n+4].i = i*g->nz2+j;
			i = round((src->s[n].x-g->x0+g->dx)/g->dx)+g->ab.np;
			src->s[n+5].i = i*g->nz2+j;
			i = round((src->s[n].x-g->x0-g->dx)/g->dx)+g->ab.np;
			j = round((src->s[n].z-g->z0+g->dz)/g->dz)+g->ab.np;
			src->s[n+6].i = i*g->nz2+j;
			i = round((src->s[n].x-g->x0      )/g->dx)+g->ab.np;
			src->s[n+7].i = i*g->nz2+j;
			i = round((src->s[n].x-g->x0+g->dx)/g->dx)+g->ab.np;
			src->s[n+8].i = i*g->nz2+j;
			
		}
	}
	fclose(fid);
}

void read_source_3d(const char filename[], const struct grid3d *g,
                    struct sourceParams *src) {
	FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    
	fscanf(fid, "%zd", &(src->nsrc));
	fscanf(fid, "%zd", &(src->nTemplate));
    if ( src->nTemplate != 1 ) {
        fprintf(stdout, "Warning, source template parameter ignored.\n");
        src->nTemplate = 1;
    }
	if ( NULL == ( src->s = (struct source *) malloc(src->nsrc*src->nTemplate*sizeof(struct source)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
	char t[5];
	for ( size_t n=0; n<src->nsrc*src->nTemplate; n+=src->nTemplate ) {
		
		fscanf(fid, "%s", t);
		if ( strcmp(t,"sx")==0 || strcmp(t,"SX")==0 || strcmp(t,"Sx")==0 ) {
			src->s[n].type = SX;
		} else if ( strcmp(t,"sz")==0 || strcmp(t,"SZ")==0 || strcmp(t,"Sz")==0 ) {
			src->s[n].type = SZ;
		} else if ( strcmp(t,"sxz")==0 || strcmp(t,"SXZ")==0 || strcmp(t,"Sxz")==0 ) {
			src->s[n].type = SXZ;
		} else if ( strcmp(t,"sf")==0 || strcmp(t,"SF")==0 || strcmp(t,"Sf")==0 ) {
			src->s[n].type = SF;
		} else if ( strcmp(t,"bulk")==0 || strcmp(t,"BULK")==0 || strcmp(t,"Bulk")==0 ) {
			src->s[n].type = BULK;
		} else {
            fprintf(stderr, "Error: type of source (%s) not implemented\n", t);
            exit(1);
        }
		
		fscanf(fid,"%lf", &(src->s[n].A) );
		fscanf(fid,"%lf", &(src->s[n].f) );
		fscanf(fid,"%lf %lf %lf", &(src->s[n].x), &(src->s[n].y), &(src->s[n].z) );
		
        //		src->s[n].A *= 1.e6;   // source strength now in Pa
		size_t i = round((src->s[n].x-g->x0)/g->dx)+g->ab.np;
		size_t j = round((src->s[n].y-g->y0)/g->dy)+g->ab.np;
        size_t k = round((src->s[n].z-g->z0)/g->dz)+g->ab.np;
		src->s[n].i = (i*g->ny2+j)*g->nz2+k;
		
	}
	fclose(fid);
}


void read_output(const char filename[], const struct grid *g, struct outputParams *out) {
	FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }

    fscanf(fid, "%zd", &(out->nrec));
	if ( NULL == ( out->r = (struct record *) malloc(out->nrec*sizeof(struct record)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }

	char t[80];
    for ( size_t n=0; n<out->nrec; ++n ) {
        
        fscanf(fid, "%s", t);
        if ( strcmp(t,"Vx")==0 || strcmp(t,"vx")==0 || strcmp(t,"VX")==0 ) {
            out->r[n].comp = VX;
        } else if ( strcmp(t,"Vy")==0 || strcmp(t,"vy")==0 || strcmp(t,"VY")==0 ) {
            out->r[n].comp = VY;
        } else if ( strcmp(t,"Vz")==0 || strcmp(t,"vz")==0 || strcmp(t,"VZ")==0 ) {
            out->r[n].comp = VZ;
        } else if ( strcmp(t,"Qx")==0 || strcmp(t,"qx")==0 || strcmp(t,"QX")==0 ) {
            out->r[n].comp = QX;
        } else if ( strcmp(t,"Qz")==0 || strcmp(t,"qz")==0 || strcmp(t,"QZ")==0 ) {
            out->r[n].comp = QZ;
        } else if ( strcmp(t,"txx")==0 || strcmp(t,"TXX")==0 ) {
            out->r[n].comp = TXX;
        } else if ( strcmp(t,"tzz")==0 || strcmp(t,"TZZ")==0 ) {
            out->r[n].comp = TZZ;
        } else if ( strcmp(t,"txy")==0 || strcmp(t,"TXY")==0 ) {
            out->r[n].comp = TXY;
        } else if ( strcmp(t,"tyz")==0 || strcmp(t,"TYZ")==0 ) {
            out->r[n].comp = TYZ;
        } else if ( strcmp(t,"txz")==0 || strcmp(t,"TXZ")==0 ) {
            out->r[n].comp = TXZ;
        } else if ( strcmp(t,"Vr")==0 || strcmp(t,"vr")==0 || strcmp(t,"VR")==0 ) {
            out->r[n].comp = VR;
        } else if ( strcmp(t,"Vtheta")==0 || strcmp(t,"vtheta")==0 || strcmp(t,"VTHETA")==0 ) {
            out->r[n].comp = VT;
        } else if ( strcmp(t,"trr")==0 || strcmp(t,"TRR")==0 ) {
            out->r[n].comp = TRR;
        } else if ( strcmp(t,"trtheta")==0 || strcmp(t,"TRTHETA")==0 ) {
            out->r[n].comp = TRT;
        } else if ( strcmp(t,"trz")==0 || strcmp(t,"TRZ")==0 ) {
            out->r[n].comp = TRZ;
        } else if ( strcmp(t,"tthetatheta")==0 || strcmp(t,"TTHETATHETA")==0 ) {
            out->r[n].comp = TTT;
        } else if ( strcmp(t,"tthetaz")==0 || strcmp(t,"TTHETAZ")==0 ) {
            out->r[n].comp = TTZ;
        } else if ( strcmp(t,"p")==0 || strcmp(t,"P")==0 ) {
            out->r[n].comp = P;
        } else if ( strcmp(t,"div")==0 || strcmp(t,"Div")==0 || strcmp(t,"DIV")==0 ) {
            out->r[n].comp = DIV;
        } else if ( strcmp(t,"curl")==0 || strcmp(t,"Curl")==0 || strcmp(t,"CURL")==0 ) {
            out->r[n].comp = CURL;
        } else {
            fprintf(stderr, "Error: particle velocity component (%s) does not exists\n", t);
            exit(1);
        }
        
        fscanf(fid, "%s", t);
        if ( strcmp(t,"Trace")==0 || strcmp(t,"trace")==0 || strcmp(t,"TRACE")==0 ) {
            out->r[n].type = TRACE;
            fscanf(fid, "%lf %lf", &(out->r[n].x), &(out->r[n].z) );
            fscanf(fid, "%lf", &(out->r[n].t0) );  // in s in the file
            fscanf(fid, "%lf", &(out->r[n].dt) );  // in ms in the file
            out->r[n].dt *= 1.e-3;
            
            size_t i = lround((out->r[n].x-g->x0)/g->dx)+g->ab.np;
            if ( g->coord == CYLINDRICAL )
                i = round((out->r[n].x-g->x0)/g->dx);

            size_t j = lround((out->r[n].z-g->z0)/g->dz)+g->ab.np;
            out->r[n].i = i*g->nz2+j;
            
        } else if ( strcmp(t,"Snapshot")==0 || strcmp(t,"snapshot")==0 || strcmp(t,"SNAPSHOT")==0 ) {
            out->r[n].type = SNAPSHOT;
            fscanf(fid, "%lf", &(out->r[n].t0) );  // in s in the file
            fscanf(fid, "%lf", &(out->r[n].dt) );  // in ms in the file
            out->r[n].dt *= 1.e-3;
            
        } else {
            fprintf(stderr, "Error: type of record (%s) not implemented\n", t);
            exit(1);
        }
        
        out->r[n].fid = NULL;
    }
    fclose(fid);
}

void read_output_3d(const char filename[], const struct grid3d *g,
                    struct outputParams *out) {
	FILE *fid;
    
    fid = fopen(filename, "r");
    if ( fid==NULL ) {
        fprintf(stderr, "Error, cannot open %s for reading\n", filename);
        exit(1);
    }
    
    fscanf(fid, "%zd", &(out->nrec));
	if ( NULL == ( out->r = (struct record *) malloc(out->nrec*sizeof(struct record)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
    
	char t[80];
    for ( size_t n=0; n<out->nrec; ++n ) {
        
        fscanf(fid, "%s", t);
        if ( strcmp(t,"Vx")==0 || strcmp(t,"vx")==0 || strcmp(t,"VX")==0 ) {
            out->r[n].comp = VX;
        } else if ( strcmp(t,"Vy")==0 || strcmp(t,"vy")==0 || strcmp(t,"VY")==0 ) {
            out->r[n].comp = VY;
        } else if ( strcmp(t,"Vz")==0 || strcmp(t,"vz")==0 || strcmp(t,"VZ")==0 ) {
            out->r[n].comp = VZ;
        } else if ( strcmp(t,"p")==0 || strcmp(t,"P")==0 ) {
            out->r[n].comp = P;
        } else {
            fprintf(stderr, "Error: particle velocity component (%s) does not exists\n", t);
            exit(1);
        }
        
        fscanf(fid, "%s", t);
        if ( strcmp(t,"Trace")==0 || strcmp(t,"trace")==0 || strcmp(t,"TRACE")==0 ) {
            out->r[n].type = TRACE;
            fscanf(fid, "%lf %lf %lf", &(out->r[n].x), &(out->r[n].y), &(out->r[n].z) );
            fscanf(fid, "%lf", &(out->r[n].t0) );  // in s in the file
            fscanf(fid, "%lf", &(out->r[n].dt) );  // in ms in the file
            out->r[n].dt *= 1.e-3;
            
            size_t i = lround((out->r[n].x-g->x0)/g->dx)+g->ab.np;
            size_t j = lround((out->r[n].y-g->y0)/g->dy)+g->ab.np;
            size_t k = lround((out->r[n].z-g->z0)/g->dz)+g->ab.np;
            out->r[n].i = (i*g->ny2+j)*g->nz2+k;
            
        } else if ( strcmp(t,"Snapshot")==0 || strcmp(t,"snapshot")==0 || strcmp(t,"SNAPSHOT")==0 ) {
            out->r[n].type = SNAPSHOT;
            fscanf(fid, "%lf", &(out->r[n].t0) );  // in s in the file
            fscanf(fid, "%lf", &(out->r[n].dt) );  // in ms in the file
            out->r[n].dt *= 1.e-3;
            
        } else {
            fprintf(stderr, "Error: type of record (%s) not implemented\n", t);
            exit(1);
        }
        
        out->r[n].fid = NULL;
    }
    fclose(fid);
}


void write_trace(const double *data, const double t, struct outputParams *out,
                 const size_t nr) {
    if ( out->r[nr].fid == NULL ) {
        struct stat sb;
        char filename[80];
        if ( stat("trc", &sb) < 0 ) {
            if ( errno == ENOENT ) {
                mkdir("trc", S_IRWXU);
            }
        }
		sprintf(filename, "trc/%s", out->basename );
        if ( stat(filename, &sb) < 0 ) {
            if ( errno == ENOENT ) {
                mkdir(filename, S_IRWXU);
            }
        }		
        sprintf(filename, "trc/%s/%s_%zd.trc", out->basename, out->basename, (nr+1) );
        out->r[nr].fid = fopen(filename, "w+");
        if ( out->r[nr].fid == NULL ) {
            fprintf(stderr, "Error, cannot open %s for output (%s)\n", filename, strerror(errno));
            exit(1);
        }
    }
    fprintf(out->r[nr].fid, "%le    %le\n", t, 1.e3*data[out->r[nr].i]);
}

void fill_data_segy(const double *data, const double t, struct outputParams *out,
                    const size_t nr) {
    size_t i = round( t/out->r[nr].dt );
    out->r[nr].data[i] = (float)data[out->r[nr].i];
}

void write_snapshot(const double *data, const double t, const struct grid *g,
                    struct outputParams *out, const size_t nr, const short plstr) {

    struct stat sb;
    char filename[80];
    if ( stat("ssh", &sb) < 0 ) {
        if ( errno == ENOENT ) {
            mkdir("ssh", S_IRWXU);
        }
    }		
	sprintf(filename, "ssh/%s", out->basename );
	if ( stat(filename, &sb) < 0 ) {
		if ( errno == ENOENT ) {
			mkdir(filename, S_IRWXU);
		}
	}		
    sprintf(filename, "ssh/%s/%s_%zd_%lf.ssh", out->basename, out->basename, (nr+1), t );
    out->r[nr].fid = fopen(filename, "w");
    if ( out->r[nr].fid == NULL ) {
        fprintf(stderr, "Error, cannot open %s for output (%s)\n", filename, strerror(errno));
        exit(1);
    }
    if ( plstr ) {
        for (size_t i=0; i<g->nx2; ++i) {
            double x = i;
            x = g->x0 + (x-g->ab.np)*g->dx;
            for (size_t j=0; j<g->nz2; ++j) {
                size_t ind = i*g->nz2 + j;
                double z = j;
                z = g->z0 + (z-g->ab.np)*g->dz;
                fprintf(out->r[nr].fid,"%le  %le  %le\n", x, z, 1.e3*data[ind]);
            }
        }
    } else {
        for (size_t i=0; i<g->nx; ++i) {
            double x = g->x0 + i*g->dx;
            for (size_t j=0; j<g->nz; ++j) {
                size_t ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
                double z = g->z0 + j*g->dz;
                fprintf(out->r[nr].fid,"%le  %le  %le\n", x, z, 1.e3*data[ind]);
            }
        }
    }
	fclose(out->r[nr].fid);
}

void write_snapshot_nc(const double *data, const double t, const struct grid *g,
                       struct outputParams *out, const size_t nr, const short plstr) {
    
    struct stat sb;
    char filename[80];
    if ( stat("ssh", &sb) < 0 ) {
        if ( errno == ENOENT ) {
            mkdir("ssh", S_IRWXU);
        }
    }
	sprintf(filename, "ssh/%s", out->basename );
	if ( stat(filename, &sb) < 0 ) {
		if ( errno == ENOENT ) {
			mkdir(filename, S_IRWXU);
		}
	}		
    sprintf(filename, "ssh/%s/%s_%zd_%lf.nc", out->basename, out->basename, (nr+1), t );
    
    int rval, ncid, z_id = -1, ids[5] = {-1,-1,-1,-1,-1}, dims[5];//, nvars, ndims;
    if (( rval=nc_create(filename, NC_CLOBBER, &ncid) ))
        ERR(rval);
    
    if ( plstr ) {
        if (( rval=nc_def_dim (ncid, "x", (size_t) g->nx2, &dims[1]) ))
            ERR(rval);
    } else {
        if (( rval=nc_def_dim (ncid, "x", (size_t) g->nx, &dims[1]) ))
            ERR(rval);
    }
    if (( rval=nc_def_var (ncid, "x", NC_DOUBLE, 1, &dims[1], &ids[1]) ))
        ERR(rval);
    
    if ( plstr ) {
        if (( rval=nc_def_dim (ncid, "y", (size_t) g->nz2, &dims[0]) ))
            ERR(rval);
    } else {
        if (( rval=nc_def_dim (ncid, "y", (size_t) g->nz, &dims[0]) ))
            ERR(rval);
    }
    if (( rval=nc_def_var (ncid, "y", NC_DOUBLE, 1, &dims[0], &ids[0]) ))
        ERR(rval);
    
    if (( rval=nc_def_var (ncid, "z", NC_FLOAT, 2, dims, &z_id) ))
        ERR(rval);

    if (( rval=nc_put_att_text (ncid, NC_GLOBAL, "Conventions", 14,  "COARDS/CF-1.0") ))
        ERR(rval);
	if (( rval=nc_put_att_text (ncid, NC_GLOBAL, "title", 8, "pve_iso") ))
		ERR(rval);
	int offset[1] = { 0 };
    if (( rval=nc_put_att_int (ncid, NC_GLOBAL, "node_offset", NC_LONG, (size_t)1, offset) ))
		ERR(rval);
	
	double dummy[2];
    if ( plstr ) {
        if ( g->coord == CARTESIAN )
            dummy[0] = g->x0 - g->ab.np*g->dx;
        else // cylindrical
            dummy[0] = g->x0;
        dummy[1] = g->x0 + (g->nx+g->ab.np-1)*g->dx;
    } else {
        dummy[0] = g->x0;
        dummy[1] = g->x0 + (g->nx-1)*g->dx;
    }
	if (( rval=nc_put_att_text (ncid, ids[1], "long_name", 2, "X") ))
        ERR(rval);
	if (( rval=nc_put_att_text (ncid, ids[1], "units", 2, "m") ))
        ERR(rval);
    if (( rval=nc_put_att_double (ncid, ids[1], "actual_range", NC_DOUBLE, (size_t)2, dummy) ))
        ERR(rval);
	
    if ( plstr ) {
        dummy[0] = g->z0 - g->ab.np*g->dz;
        dummy[1] = g->z0 + (g->nz+g->ab.np-1)*g->dz;
    } else {
        dummy[0] = g->z0;
        dummy[1] = g->z0 + (g->nz-1)*g->dz;
    }
	if (( rval=nc_put_att_text (ncid, ids[0], "long_name", 2, "Y") ))
        ERR(rval);
	if (( rval=nc_put_att_text (ncid, ids[0], "units", 2, "m") ))
        ERR(rval);
    if (( rval=nc_put_att_double (ncid, ids[0], "actual_range", NC_DOUBLE, (size_t)2, dummy) ))
        ERR(rval);
	
	double scale[1] = { 1. };
	if (( rval=nc_put_att_text (ncid, z_id, "units", 5, "mm/s") ))
        ERR(rval);
	if (( rval=nc_put_att_double (ncid, z_id, "scale_factor", NC_DOUBLE, (size_t)1, scale) ))
        ERR(rval);
	
	double nan[1] = { fabs( 0./0. ) };
	if (( rval=nc_put_att_double (ncid, z_id, "_FillValue", NC_FLOAT, (size_t) 1, nan) ))
        ERR(rval);
	
    dummy[0] = dummy[1] = 1.e3*data[0];
    if ( plstr ) {
        for ( size_t ind=1; ind<g->nx2*g->nz2; ++ind ) {
            dummy[0] = dummy[0] < 1.e3*data[ind] ? dummy[0] : 1.e3*data[ind];
            dummy[1] = dummy[1] > 1.e3*data[ind] ? dummy[1] : 1.e3*data[ind];
        }
    } else {
        for (size_t j=0; j<g->nz; ++j) {
            for (size_t i=0; i<g->nx; ++i) {
                size_t ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
                dummy[0] = dummy[0] < 1.e3*data[ind] ? dummy[0] : 1.e3*data[ind];
                dummy[1] = dummy[1] > 1.e3*data[ind] ? dummy[1] : 1.e3*data[ind];
            }
        }
    }
	if (( rval=nc_put_att_double (ncid, z_id, "actual_range", NC_DOUBLE, (size_t)2, dummy) ))
        ERR(rval);
	
	if (( rval=nc_enddef (ncid) ))
        ERR(rval);
	
	double *x, *y;
	float *tmp_f;
	size_t start[2] = {0,0}, edge[2] = {1,1};
    if ( plstr ) {
        if ( NULL == ( tmp_f = (float*) malloc( g->nx2*sizeof(float) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( x = (double*) malloc( g->nx2*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( y = (double*) malloc( g->nz2*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        for (size_t i=0; i<g->nx2; ++i) {
            if ( g->coord == CARTESIAN )
                x[i] = g->x0 + ((1.0*i)-g->ab.np)*g->dx;
            else // cylindrical
                x[i] = g->x0 + i*g->dx;
        }
        if (( rval=nc_put_var_double (ncid, ids[1], x) ))
            ERR(rval);
        for (size_t i=0; i<g->nz2; ++i) {
            y[i] = g->z0 + ((1.0*i)-g->ab.np)*g->dz;
        }
        if (( rval=nc_put_var_double (ncid, ids[0], y) ))
            ERR(rval);
        
        edge[1] = g->nx2;
        for (size_t j=0; j<g->nz2; ++j) {
            start[0] = j;
            for (size_t i=0; i<g->nx2; ++i) {
                size_t ind = i*g->nz2 + j;
                tmp_f[i] = (float)data[ind]*1.e3f;
            }
            if (( rval=nc_put_vara_float (ncid, z_id, start, edge, tmp_f) ))
                ERR(rval);
        }
    } else {
        if ( NULL == ( tmp_f = (float*) malloc( g->nx*sizeof(float) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( x = (double*) malloc( g->nx*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( y = (double*) malloc( g->nz*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        for (size_t i=0; i<g->nx; ++i) {
            x[i] = g->x0 + i*g->dx;
        }
        if (( rval=nc_put_var_double (ncid, ids[1], x) ))
            ERR(rval);
        for (size_t i=0; i<g->nz; ++i) {
            y[i] = g->z0 + i*g->dz;
        }
        if (( rval=nc_put_var_double (ncid, ids[0], y) ))
            ERR(rval);
        
        edge[1] = g->nx;
        for (size_t j=0; j<g->nz; ++j) {
            start[0] = j;
            for (size_t i=0; i<g->nx; ++i) {
                size_t ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
                tmp_f[i] = (float)data[ind]*1.e3f;
            }
            if (( rval=nc_put_vara_float (ncid, z_id, start, edge, tmp_f) ))
                ERR(rval);
        }        
    }

	if (( rval=nc_close (ncid) ))
        ERR(rval);

	free(x);
	free(y);
	free(tmp_f);	
}


void write_snapshot3d_nc(const double *data, const double t,
                         const struct grid3d *g, struct outputParams *out,
                         const size_t nr, const short plstr) {
    
    struct stat sb;
    if ( stat("ssh", &sb) < 0 ) {
        if ( errno == ENOENT ) {
            mkdir("ssh", S_IRWXU);
        }
    }
    char filename[80];
	sprintf(filename, "ssh/%s", out->basename );
	if ( stat(filename, &sb) < 0 ) {
		if ( errno == ENOENT ) {
			mkdir(filename, S_IRWXU);
		}
	}
    sprintf(filename, "ssh/%s/%s_%zd_%lf.nc", out->basename, out->basename, (nr+1), t );
    
    int rval, ncid, z_id = -1, ids[5] = {-1,-1,-1,-1,-1}, dims[5];//, nvars, ndims;
    if (( rval=nc_create(filename, NC_CLOBBER, &ncid) ))
        ERR(rval);
    
    if ( plstr ) {
        if (( rval=nc_def_dim (ncid, "x", (size_t) g->nx2, &dims[1]) ))
            ERR(rval);
    } else {
        if (( rval=nc_def_dim (ncid, "x", (size_t) g->nx, &dims[1]) ))
            ERR(rval);
    }
    if (( rval=nc_def_var (ncid, "x", NC_DOUBLE, 1, &dims[1], &ids[1]) ))
        ERR(rval);
    
    if ( plstr ) {
        if (( rval=nc_def_dim (ncid, "y", (size_t) g->nz2, &dims[0]) ))
            ERR(rval);
    } else {
        if (( rval=nc_def_dim (ncid, "y", (size_t) g->nz, &dims[0]) ))
            ERR(rval);
    }
    if (( rval=nc_def_var (ncid, "y", NC_DOUBLE, 1, &dims[0], &ids[0]) ))
        ERR(rval);
    
    if (( rval=nc_def_var (ncid, "z", NC_FLOAT, 2, dims, &z_id) ))
        ERR(rval);
    
    if (( rval=nc_put_att_text (ncid, NC_GLOBAL, "Conventions", 14,  "COARDS/CF-1.0") ))
        ERR(rval);
	if (( rval=nc_put_att_text (ncid, NC_GLOBAL, "title", 8, "pve_iso") ))
		ERR(rval);
	int offset[1] = { 0 };
    if (( rval=nc_put_att_int (ncid, NC_GLOBAL, "node_offset", NC_LONG, (size_t)1, offset) ))
		ERR(rval);
	
	double dummy[2];
    if ( plstr ) {
        dummy[0] = g->x0 - g->ab.np*g->dx;
        dummy[1] = g->x0 + (g->nx+g->ab.np-1)*g->dx;
    } else {
        dummy[0] = g->x0;
        dummy[1] = g->x0 + (g->nx-1)*g->dx;
    }
	if (( rval=nc_put_att_text (ncid, ids[1], "long_name", 2, "X") ))
        ERR(rval);
	if (( rval=nc_put_att_text (ncid, ids[1], "units", 2, "m") ))
        ERR(rval);
    if (( rval=nc_put_att_double (ncid, ids[1], "actual_range", NC_DOUBLE, (size_t)2, dummy) ))
        ERR(rval);
	
    if ( plstr ) {
        dummy[0] = g->z0 - g->ab.np*g->dz;
        dummy[1] = g->z0 + (g->nz+g->ab.np-1)*g->dz;
    } else {
        dummy[0] = g->z0;
        dummy[1] = g->z0 + (g->nz-1)*g->dz;
    }
	if (( rval=nc_put_att_text (ncid, ids[0], "long_name", 2, "Y") ))
        ERR(rval);
	if (( rval=nc_put_att_text (ncid, ids[0], "units", 2, "m") ))
        ERR(rval);
    if (( rval=nc_put_att_double (ncid, ids[0], "actual_range", NC_DOUBLE, (size_t)2, dummy) ))
        ERR(rval);
	
	double scale[1] = { 1. };
	if (( rval=nc_put_att_text (ncid, z_id, "units", 5, "mm/s") ))
        ERR(rval);
	if (( rval=nc_put_att_double (ncid, z_id, "scale_factor", NC_DOUBLE, (size_t)1, scale) ))
        ERR(rval);
	
	double nan[1] = { fabs( 0./0. ) };
	if (( rval=nc_put_att_double (ncid, z_id, "_FillValue", NC_FLOAT, (size_t) 1, nan) ))
        ERR(rval);
	
    dummy[0] = dummy[1] = 1.e3*data[0];
    if ( plstr ) {
        for ( size_t ind=1; ind<g->nx2*g->nz2; ++ind ) {
            dummy[0] = dummy[0] < 1.e3*data[ind] ? dummy[0] : 1.e3*data[ind];
            dummy[1] = dummy[1] > 1.e3*data[ind] ? dummy[1] : 1.e3*data[ind];
        }
    } else {
        for (size_t j=0; j<g->nz; ++j) {
            for (size_t i=0; i<g->nx; ++i) {
                size_t ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
                dummy[0] = dummy[0] < 1.e3*data[ind] ? dummy[0] : 1.e3*data[ind];
                dummy[1] = dummy[1] > 1.e3*data[ind] ? dummy[1] : 1.e3*data[ind];
            }
        }
    }
	if (( rval=nc_put_att_double (ncid, z_id, "actual_range", NC_DOUBLE, (size_t)2, dummy) ))
        ERR(rval);
	
	if (( rval=nc_enddef (ncid) ))
        ERR(rval);
	
	double *x, *y;
	float *tmp_f;
	size_t start[2] = {0,0}, edge[2] = {1,1};
    if ( plstr ) {
        if ( NULL == ( tmp_f = (float*) malloc( g->nx2*sizeof(float) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( x = (double*) malloc( g->nx2*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( y = (double*) malloc( g->nz2*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        for (size_t i=0; i<g->nx2; ++i) {
            x[i] = g->x0 + ((1.0*i)-g->ab.np)*g->dx;
        }
        if (( rval=nc_put_var_double (ncid, ids[1], x) ))
            ERR(rval);
        for (size_t i=0; i<g->nz2; ++i) {
            y[i] = g->z0 + ((1.0*i)-g->ab.np)*g->dz;
        }
        if (( rval=nc_put_var_double (ncid, ids[0], y) ))
            ERR(rval);
        
        edge[1] = g->nx2;
        for (size_t j=0; j<g->nz2; ++j) {
            start[0] = j;
            for (size_t i=0; i<g->nx2; ++i) {
                size_t ind = i*g->nz2 + j;
                tmp_f[i] = (float)data[ind]*1.e3f;
            }
            if (( rval=nc_put_vara_float (ncid, z_id, start, edge, tmp_f) ))
                ERR(rval);
        }
    } else {
        if ( NULL == ( tmp_f = (float*) malloc( g->nx*sizeof(float) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( x = (double*) malloc( g->nx*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( y = (double*) malloc( g->nz*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        for (size_t i=0; i<g->nx; ++i) {
            x[i] = g->x0 + i*g->dx;
        }
        if (( rval=nc_put_var_double (ncid, ids[1], x) ))
            ERR(rval);
        for (size_t i=0; i<g->nz; ++i) {
            y[i] = g->z0 + i*g->dz;
        }
        if (( rval=nc_put_var_double (ncid, ids[0], y) ))
            ERR(rval);
        
        edge[1] = g->nx;
        for (size_t j=0; j<g->nz; ++j) {
            start[0] = j;
            for (size_t i=0; i<g->nx; ++i) {
                size_t ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
                tmp_f[i] = (float)data[ind]*1.e3f;
            }
            if (( rval=nc_put_vara_float (ncid, z_id, start, edge, tmp_f) ))
                ERR(rval);
        }
    }
    
	if (( rval=nc_close (ncid) ))
        ERR(rval);
    
	free(x);
	free(y);
	free(tmp_f);
}


void write_field_nc(const double *data, const char *fieldname,
                    const char *units, const struct grid *g,
                    struct outputParams *out, const short plstr) {
    char filename[80];
    sprintf(filename, "%s_%s.nc", out->basename, fieldname );
    
    int rval, ncid, z_id = -1, ids[5] = {-1,-1,-1,-1,-1}, dims[5];//, nvars, ndims;
    if (( rval=nc_create(filename, NC_CLOBBER, &ncid) ))
        ERR(rval);
    
    if ( plstr ) {
        if (( rval=nc_def_dim (ncid, "x", (size_t) g->nx2, &dims[1]) ))
            ERR(rval);
    } else {
        if (( rval=nc_def_dim (ncid, "x", (size_t) g->nx, &dims[1]) ))
            ERR(rval);
    }
    if (( rval=nc_def_var (ncid, "x", NC_DOUBLE, 1, &dims[1], &ids[1]) ))
        ERR(rval);
    
    if ( plstr ) {
        if (( rval=nc_def_dim (ncid, "y", (size_t) g->nz2, &dims[0]) ))
            ERR(rval);
    } else {
        if (( rval=nc_def_dim (ncid, "y", (size_t) g->nz, &dims[0]) ))
            ERR(rval);
    }
    if (( rval=nc_def_var (ncid, "y", NC_DOUBLE, 1, &dims[0], &ids[0]) ))
        ERR(rval);
    
    if (( rval=nc_def_var (ncid, "z", NC_FLOAT, 2, dims, &z_id) ))
        ERR(rval);
    
    if (( rval=nc_put_att_text (ncid, NC_GLOBAL, "Conventions", 14,  "COARDS/CF-1.0") ))
        ERR(rval);
	if (( rval=nc_put_att_text (ncid, NC_GLOBAL, "title", 8, "pve_iso") ))
		ERR(rval);
	int offset[1] = { 0 };
    if (( rval=nc_put_att_int (ncid, NC_GLOBAL, "node_offset", NC_LONG, (size_t)1, offset) ))
		ERR(rval);
	
	double dummy[2];
    if ( plstr ) {
        if ( g->coord == CARTESIAN )
            dummy[0] = g->x0 - g->ab.np*g->dx;
        else // cylindrical
            dummy[0] = g->x0;
        dummy[1] = g->x0 + (g->nx+g->ab.np-1)*g->dx;
    } else {
        dummy[0] = g->x0;
        dummy[1] = g->x0 + (g->nx-1)*g->dx;
    }
	if (( rval=nc_put_att_text (ncid, ids[1], "long_name", 2, "X") ))
        ERR(rval);
	if (( rval=nc_put_att_text (ncid, ids[1], "units", 2, "m") ))
        ERR(rval);
    if (( rval=nc_put_att_double (ncid, ids[1], "actual_range", NC_DOUBLE, (size_t)2, dummy) ))
        ERR(rval);
	
    if ( plstr ) {
        dummy[0] = g->z0 - g->ab.np*g->dz;
        dummy[1] = g->z0 + (g->nz+g->ab.np-1)*g->dz;
    } else {
        dummy[0] = g->z0;
        dummy[1] = g->z0 + (g->nz-1)*g->dz;
    }
	if (( rval=nc_put_att_text (ncid, ids[0], "long_name", 2, "Y") ))
        ERR(rval);
	if (( rval=nc_put_att_text (ncid, ids[0], "units", 2, "m") ))
        ERR(rval);
    if (( rval=nc_put_att_double (ncid, ids[0], "actual_range", NC_DOUBLE, (size_t)2, dummy) ))
        ERR(rval);
	
	double scale[1] = { 1. };
	if (( rval=nc_put_att_text (ncid, z_id, "units", strlen(units), units) ))
        ERR(rval);
	if (( rval=nc_put_att_double (ncid, z_id, "scale_factor", NC_DOUBLE, (size_t)1, scale) ))
        ERR(rval);
	
	double nan[1] = { fabs( 0./0. ) };
	if (( rval=nc_put_att_double (ncid, z_id, "_FillValue", NC_FLOAT, (size_t) 1, nan) ))
        ERR(rval);
	
    dummy[0] = dummy[1] = data[0];
    if ( plstr ) {
        for ( size_t ind=1; ind<g->nx2*g->nz2; ++ind ) {
            dummy[0] = dummy[0] < data[ind] ? dummy[0] : data[ind];
            dummy[1] = dummy[1] > data[ind] ? dummy[1] : data[ind];
        }
    } else {
        for (size_t j=0; j<g->nz; ++j) {
            for (size_t i=0; i<g->nx; ++i) {
                size_t ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
                dummy[0] = dummy[0] < data[ind] ? dummy[0] : data[ind];
                dummy[1] = dummy[1] > data[ind] ? dummy[1] : data[ind];
            }
        }
    }
	if (( rval=nc_put_att_double (ncid, z_id, "actual_range", NC_DOUBLE, (size_t)2, dummy) ))
        ERR(rval);
	
	if (( rval=nc_enddef (ncid) ))
        ERR(rval);
	
	double *x, *y;
	float *tmp_f;
	size_t start[2] = {0,0}, edge[2] = {1,1};
    if ( plstr ) {
        if ( NULL == ( tmp_f = (float*) malloc( g->nx2*sizeof(float) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( x = (double*) malloc( g->nx2*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( y = (double*) malloc( g->nz2*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        for (size_t i=0; i<g->nx2; ++i) {
            if ( g->coord == CARTESIAN )
                x[i] = g->x0 + ((1.0*i)-g->ab.np)*g->dx;
            else // cylindrical
                x[i] = g->x0 + i*g->dx;
        }
        if (( rval=nc_put_var_double (ncid, ids[1], x) ))
            ERR(rval);
        for (size_t i=0; i<g->nz2; ++i) {
            y[i] = g->z0 + ((1.0*i)-g->ab.np)*g->dz;
        }
        if (( rval=nc_put_var_double (ncid, ids[0], y) ))
            ERR(rval);
        
        edge[1] = g->nx2;
        for (size_t j=0; j<g->nz2; ++j) {
            start[0] = j;
            for (size_t i=0; i<g->nx2; ++i) {
                size_t ind = i*g->nz2 + j;
                tmp_f[i] = (float)data[ind];
            }
            if (( rval=nc_put_vara_float (ncid, z_id, start, edge, tmp_f) ))
                ERR(rval);
        }
    } else {
        if ( NULL == ( tmp_f = (float*) malloc( g->nx*sizeof(float) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( x = (double*) malloc( g->nx*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        if ( NULL == ( y = (double*) malloc( g->nz*sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
        for (size_t i=0; i<g->nx; ++i) {
            x[i] = g->x0 + i*g->dx;
        }
        if (( rval=nc_put_var_double (ncid, ids[1], x) ))
            ERR(rval);
        for (size_t i=0; i<g->nz; ++i) {
            y[i] = g->z0 + i*g->dz;
        }
        if (( rval=nc_put_var_double (ncid, ids[0], y) ))
            ERR(rval);
        
        edge[1] = g->nx;
        for (size_t j=0; j<g->nz; ++j) {
            start[0] = j;
            for (size_t i=0; i<g->nx; ++i) {
                size_t ind = (i+g->ab.np)*g->nz2 + j+g->ab.np;
                tmp_f[i] = (float)data[ind];
            }
            if (( rval=nc_put_vara_float (ncid, z_id, start, edge, tmp_f) ))
                ERR(rval);
        }
    }
    
	if (( rval=nc_close (ncid) ))
        ERR(rval);
    
	free(x);
	free(y);
	free(tmp_f);
}



void save_checkpoint(const size_t it, const struct grid *g,
					 const struct inputParams *p,
					 const struct sourceParams *s,
					 const struct outputParams *out,
					 const struct saveEnergy *se,
					 const struct fac_pml *fp1,
					 const struct fac_pml *fp23,
					 const struct fac_pml *fp4,
					 const struct mem_pml mem[3],
					 const struct computationVariables *c,
					 const double *coeff_v_1,
					 const double *coeff_v_3,    
					 const double *coeff_q_1,
					 const double *coeff_q_3,
					 const double *W,
					 const double *Ws,
					 const double *Wtmp,
					 const double *delta		 
					 ) {
	
	struct stat sb;
    if ( stat("chkpt", &sb) < 0 ) {
        if ( errno == ENOENT ) {
            mkdir("chkpt", S_IRWXU);
        }
    }
    char filename[80];
    sprintf(filename, "chkpt/%s_%zd.dat", p->basename, it );
	FILE *fid = fopen(filename, "w");
    if ( fid == NULL ) {
        fprintf(stderr, "Error, cannot open %s for output (%s)\n", filename, strerror(errno));
        exit(1);
    }
	
	fwrite(&it, sizeof(size_t), 1, fid);
	fwrite(g, sizeof(struct grid), 1, fid);
	fwrite(p, sizeof(struct inputParams), 1, fid);
	fwrite(&(s->nsrc), sizeof(size_t), 1, fid);
	fwrite(&(s->nTemplate), sizeof(size_t), 1, fid);
	for ( size_t n=0; n<s->nsrc*s->nTemplate; ++n ) {
		fwrite(&(s->s[n]), sizeof(struct source), 1, fid);
		fwrite(s->s[n].fct, sizeof(double), s->s[n].length, fid);
	}
	fwrite(out, sizeof(struct outputParams), 1, fid);
	for ( size_t nr=0; nr<out->nrec; ++nr ) {
		fwrite(&(out->r[nr]), sizeof(struct record), 1, fid);
        
        if ( out->r[nr].type == TRACE && p->segy==0 ) {
            size_t nsamples = round( p->duration / out->r[nr].dt );
            fwrite(out->r[nr].data, sizeof(float), nsamples, fid);
		}
        else if ( out->r[nr].type == TRACE && p->segy==0 ) {
			fflush(out->r[nr].fid);
			rewind(out->r[nr].fid);
			
			if ( stat("chkpt/trc", &sb) < 0 ) {
				if ( errno == ENOENT ) {
					mkdir("chkpt/trc", S_IRWXU);
				}
			}
			char filename2[100];
			sprintf(filename2, "chkpt/trc/%s_%zd", out->basename, it );
			if ( stat(filename2, &sb) < 0 ) {
				if ( errno == ENOENT ) {
					mkdir(filename2, S_IRWXU);
				}
			}		
			sprintf(filename2, "chkpt/trc/%s_%zd/%s_%zd.trc", out->basename, it, out->basename, (nr+1) );
			FILE *fid2 = fopen(filename2, "w");
			if ( fid2 == NULL ) {
				fprintf(stderr, "Error, cannot open %s for output (%s)\n", filename2, strerror(errno));
				exit(1);
			}
			char ch;
			while(!feof(out->r[nr].fid)) {
				ch = fgetc(out->r[nr].fid);
				if(ferror(out->r[nr].fid)) {
					fprintf(stderr, "Error reading source file (trace).\n");
					exit(1);
				}
				if(!feof(out->r[nr].fid)) fputc(ch, fid2);
				if(ferror(fid2)) {
					fprintf(stderr, "Error writing destination file (trace).\n");
					exit(1);
				}
			}
			fclose(fid2);
		}
	}
	if ( p->saveEnergy ) {
		fwrite(se, sizeof(struct saveEnergy), 1, fid);
		size_t npts = (se->i2E-se->i1E+1)*(se->j2E-se->j1E+1);
		fwrite(se->rho11, sizeof(double), npts, fid);
		fwrite(se->rho12, sizeof(double), npts, fid);
		fwrite(se->rho22, sizeof(double), npts, fid);
		
		struct stat sb;
		if ( stat("chkpt/E", &sb) < 0 ) {
			if ( errno == ENOENT ) {
				mkdir("chkpt/E", S_IRWXU);
			}
		}
		char fname[100];
		sprintf(fname, "chkpt/E/%s_%zd_E.dat", p->basename, it );
		FILE *fid2 = fopen(fname, "w");
		if ( fid2 == NULL ) {
			fprintf(stderr, "Error, cannot open %s for output (%s)\n", fname, strerror(errno));
			exit(1);
		}
		char ch;
		fflush(se->fid);
		rewind(se->fid);
		while(!feof(se->fid)) {
			ch = fgetc(se->fid);
			if(ferror(se->fid)) {
				fprintf(stderr, "Error reading source file (energy).\n");
				exit(1);
			}
			if(!feof(se->fid)) fputc(ch, fid2);
			if(ferror(fid2)) {
				fprintf(stderr, "Error writing destination file (energy).\n");
				exit(1);
			}
		}
		fclose(fid2);
	}
	
	fwrite(fp1->k_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->k_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->kh_x, sizeof(double), g->ab.np, fid);
	fwrite(fp1->kh_z, sizeof(double), g->ab.np, fid);
    fwrite(fp1->kH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp1->kH_z, sizeof(double), g->ab.np+1, fid);
    fwrite(fp1->b_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->c_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->bh_x, sizeof(double), g->ab.np, fid);
    fwrite(fp1->ch_x, sizeof(double), g->ab.np, fid);
    fwrite(fp1->bH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp1->cH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp1->b_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->c_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->bh_z, sizeof(double), g->ab.np, fid);
    fwrite(fp1->ch_z, sizeof(double), g->ab.np, fid);
    fwrite(fp1->bH_z, sizeof(double), g->ab.np+1, fid);
    fwrite(fp1->cH_z, sizeof(double), g->ab.np+1, fid);
	
    fwrite(fp23->b_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp23->c_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp23->bh_x, sizeof(double), g->ab.np, fid);
    fwrite(fp23->ch_x, sizeof(double), g->ab.np, fid);
    fwrite(fp23->bH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp23->cH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp23->b_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp23->c_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp23->bh_z, sizeof(double), g->ab.np, fid);
    fwrite(fp23->ch_z, sizeof(double), g->ab.np, fid);
    fwrite(fp23->bH_z, sizeof(double), g->ab.np+1, fid);
    fwrite(fp23->cH_z, sizeof(double), g->ab.np+1, fid);
	
    fwrite(fp4->b_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp4->c_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp4->bh_x, sizeof(double), g->ab.np, fid);
    fwrite(fp4->ch_x, sizeof(double), g->ab.np, fid);
    fwrite(fp4->bH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp4->cH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp4->b_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp4->c_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp4->bh_z, sizeof(double), g->ab.np, fid);
    fwrite(fp4->ch_z, sizeof(double), g->ab.np, fid);
    fwrite(fp4->bH_z, sizeof(double), g->ab.np+1, fid);
    fwrite(fp4->cH_z, sizeof(double), g->ab.np+1, fid);
	
	for ( size_t n=0; n<3; ++n ) {
		fwrite(mem[n].dx_txx, sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dx_p  , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dx_txz, sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dx_vx , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dx_qx , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dx_vz , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dz_txz, sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fwrite(mem[n].dz_tzz, sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fwrite(mem[n].dz_p  , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fwrite(mem[n].dz_vz , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fwrite(mem[n].dz_qz , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fwrite(mem[n].dz_vx , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);		
	}
	
	size_t nnodes = g->nx2 * g->nz2;
    fwrite(c->mu     , sizeof(double), nnodes, fid);
    fwrite(c->epsilon, sizeof(double), nnodes, fid);
    fwrite(c->E      , sizeof(double), nnodes, fid);
    fwrite(c->M      , sizeof(double), nnodes, fid);
    fwrite(c->alpha  , sizeof(double), nnodes, fid);
    fwrite(c->varphi , sizeof(double), nnodes, fid);
    fwrite(c->tau_s  , sizeof(double), nnodes, fid);
    fwrite(c->rho_i  , sizeof(double), nnodes, fid);
    fwrite(c->rho_j  , sizeof(double), nnodes, fid);
    fwrite(c->rho_f_i, sizeof(double), nnodes, fid);
    fwrite(c->rho_f_j, sizeof(double), nnodes, fid);
    fwrite(c->nk_i   , sizeof(double), nnodes, fid);
    fwrite(c->nk_j   , sizeof(double), nnodes, fid);
    fwrite(c->mu_ij  , sizeof(double), nnodes, fid);
    fwrite(c->m_i    , sizeof(double), nnodes, fid);
    fwrite(c->m_j    , sizeof(double), nnodes, fid);
	
	fwrite(coeff_v_1, sizeof(double), nnodes, fid);
	fwrite(coeff_v_3, sizeof(double), nnodes, fid);
	fwrite(coeff_q_1, sizeof(double), nnodes, fid);
	fwrite(coeff_q_3, sizeof(double), nnodes, fid);
	fwrite(W, sizeof(double), 9*nnodes, fid);
	fwrite(Ws, sizeof(double), 9*nnodes, fid);
	fwrite(Wtmp, sizeof(double), 9*nnodes, fid);
	fwrite(delta, sizeof(double), 9*nnodes, fid);
    
	fclose(fid);
}

void read_checkpoint1(const char *filename, size_t *it, struct grid *g,
					  struct inputParams *p) {
	char   checkpointfile[80];
	strcpy(checkpointfile, filename);
	FILE *fid = fopen(filename, "r");
    if ( fid == NULL ) {
        fprintf(stderr, "Error, cannot open %s for input (%s)\n", filename, strerror(errno));
        exit(1);
    }
	fread(it, sizeof(size_t), 1, fid);
	fread(g, sizeof(struct grid), 1, fid);
	fread(p, sizeof(struct inputParams), 1, fid);
	strcpy(p->checkpointfile, checkpointfile);
	fclose(fid);
}

void read_checkpoint2(const char *filename, const struct grid *g,
						 const struct inputParams *p, struct sourceParams *s,
						 struct outputParams *out, struct saveEnergy *se,
						 struct fac_pml *fp1, struct fac_pml *fp23,
						 struct fac_pml *fp4, struct mem_pml mem[3],
						 struct computationVariables *c, double *coeff_v_i,
						 double *coeff_v_j, double *coeff_q_i, double *coeff_q_j,
						 double *W, double *Ws, double *Wtmp, double *delta) {
	
	FILE *fid = fopen(filename, "r");
    if ( fid == NULL ) {
        fprintf(stderr, "Error, cannot open %s for input (%s)\n", filename, strerror(errno));
        exit(1);
    }
	size_t it;
	fread(&it, sizeof(size_t), 1, fid);
	
	fseek(fid, sizeof(struct grid)+sizeof(struct inputParams), SEEK_CUR);
	
	fread(&(s->nsrc), sizeof(size_t), 1, fid);
	fread(&(s->nTemplate), sizeof(size_t), 1, fid);
	if ( NULL == ( s->s = (struct source *) malloc(s->nsrc*s->nTemplate*sizeof(struct source)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	for ( size_t n=0; n<s->nsrc*s->nTemplate; ++n ) {
		fread(&(s->s[n]), sizeof(struct source), 1, fid);
		if ( NULL == ( s->s[n].fct = malloc( s->s[n].length * sizeof(double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		fread(s->s[n].fct, sizeof(double), s->s[n].length, fid);
	}
	fread(out, sizeof(struct outputParams), 1, fid);
	if ( NULL == ( out->r = (struct record *) malloc(out->nrec*sizeof(struct record)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	for ( size_t nr=0; nr<out->nrec; ++nr ) {
		fread(&(out->r[nr]), sizeof(struct record), 1, fid);
		
        if ( out->r[nr].type == TRACE && p->segy == 1 ) {
            size_t nsamples = round( p->duration / out->r[nr].dt );
            if ( NULL == ( out->r[nr].data = (float *) malloc( nsamples * sizeof(float) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            fread(out->r[nr].data, sizeof(float), nsamples, fid);
        }
        else if ( out->r[nr].type == TRACE && p->segy == 0 ) {
			
            struct stat sb;
            if ( stat("trc", &sb) < 0 ) {
                if ( errno == ENOENT ) {
                    mkdir("trc", S_IRWXU);
                }
            }
            char filename1[100];
            char filename2[80];
            sprintf(filename2, "trc/%s", out->basename );
            if ( stat(filename2, &sb) < 0 ) {
                if ( errno == ENOENT ) {
                    mkdir(filename2, S_IRWXU);
                }
            }
            sprintf(filename1, "chkpt/trc/%s_%zd/%s_%zd.trc", out->basename, it, out->basename, (nr+1) );
            sprintf(filename2, "trc/%s/%s_%zd.trc", out->basename, out->basename, (nr+1) );
            FILE *fid1 = fopen(filename1, "r");
            if ( fid1 == NULL ) {
                fprintf(stderr, "Error, cannot open %s (%s)\n", filename1, strerror(errno));
                exit(1);
            }
            out->r[nr].fid = fopen(filename2, "w+");
            if ( out->r[nr].fid == NULL ) {
                fprintf(stderr, "Error, cannot open %s (%s)\n", filename2, strerror(errno));
                exit(1);
            }
            char ch;
            while(!feof(fid1)) {
                ch = fgetc(fid1);
                if(ferror(fid1)) {
                    fprintf(stderr, "Error reading source file (trace).\n");
                    exit(1);
                }
                if(!feof(fid1)) fputc(ch, out->r[nr].fid);
                if(ferror(out->r[nr].fid)) {
                    fprintf(stderr, "Error writing destination file (trace).\n");
                    exit(1);
                }
            }
            fclose(fid1);
        }
	}
	if ( p->saveEnergy ) {
		fread(se, sizeof(struct saveEnergy), 1, fid);
		size_t npts = (se->i2E-se->i1E+1)*(se->j2E-se->j1E+1);
		if ( NULL == ( se->rho11 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( se->rho12 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( se->rho22 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		fread(se->rho11, sizeof(double), npts, fid);
		fread(se->rho12, sizeof(double), npts, fid);
		fread(se->rho22, sizeof(double), npts, fid);
		
		struct stat sb;
		if ( stat("E", &sb) < 0 ) {
			if ( errno == ENOENT ) {
				mkdir("E", S_IRWXU);
			}
		}
		char fname1[100];
		char fname[100];
		sprintf(fname1, "chkpt/E/%s_%zd_E.dat", p->basename, it );
		sprintf(fname, "E/%s_E.dat", p->basename );
		FILE *fid1 = fopen(fname1, "r");
		if ( fid1 == NULL ) {
			fprintf(stderr, "Error, cannot open %s (%s)\n", fname1, strerror(errno));
			exit(1);
		}
		se->fid = fopen(fname, "w+");
		if ( se->fid == NULL ) {
			fprintf(stderr, "Error, cannot open %s for output (%s)\n", fname, strerror(errno));
			exit(1);
		}
		char ch;
		while(!feof(fid1)) {
			ch = fgetc(fid1);
			if(ferror(fid1)) {
				fprintf(stderr, "Error reading source file (energy).\n");
				exit(1);
			}
			if(!feof(fid1)) fputc(ch, se->fid);
			if(ferror(se->fid)) {
				fprintf(stderr, "Error writing destination file (energy).\n");
				exit(1);
			}
		}
		fclose(fid1);
	}
	
	fread(fp1->k_x,  sizeof(double), g->ab.np, fid);
    fread(fp1->k_z,  sizeof(double), g->ab.np, fid);
    fread(fp1->kh_x, sizeof(double), g->ab.np, fid);
	fread(fp1->kh_z, sizeof(double), g->ab.np, fid);
    fread(fp1->kH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp1->kH_z, sizeof(double), g->ab.np+1, fid);
    fread(fp1->b_x,  sizeof(double), g->ab.np, fid);
    fread(fp1->c_x,  sizeof(double), g->ab.np, fid);
    fread(fp1->bh_x, sizeof(double), g->ab.np, fid);
    fread(fp1->ch_x, sizeof(double), g->ab.np, fid);
    fread(fp1->bH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp1->cH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp1->b_z,  sizeof(double), g->ab.np, fid);
    fread(fp1->c_z,  sizeof(double), g->ab.np, fid);
    fread(fp1->bh_z, sizeof(double), g->ab.np, fid);
    fread(fp1->ch_z, sizeof(double), g->ab.np, fid);
    fread(fp1->bH_z, sizeof(double), g->ab.np+1, fid);
    fread(fp1->cH_z, sizeof(double), g->ab.np+1, fid);
	
    fread(fp23->b_x,  sizeof(double), g->ab.np, fid);
    fread(fp23->c_x,  sizeof(double), g->ab.np, fid);
    fread(fp23->bh_x, sizeof(double), g->ab.np, fid);
    fread(fp23->ch_x, sizeof(double), g->ab.np, fid);
    fread(fp23->bH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp23->cH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp23->b_z,  sizeof(double), g->ab.np, fid);
    fread(fp23->c_z,  sizeof(double), g->ab.np, fid);
    fread(fp23->bh_z, sizeof(double), g->ab.np, fid);
    fread(fp23->ch_z, sizeof(double), g->ab.np, fid);
    fread(fp23->bH_z, sizeof(double), g->ab.np+1, fid);
    fread(fp23->cH_z, sizeof(double), g->ab.np+1, fid);
	
    fread(fp4->b_x,  sizeof(double), g->ab.np, fid);
    fread(fp4->c_x,  sizeof(double), g->ab.np, fid);
    fread(fp4->bh_x, sizeof(double), g->ab.np, fid);
    fread(fp4->ch_x, sizeof(double), g->ab.np, fid);
    fread(fp4->bH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp4->cH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp4->b_z,  sizeof(double), g->ab.np, fid);
    fread(fp4->c_z,  sizeof(double), g->ab.np, fid);
    fread(fp4->bh_z, sizeof(double), g->ab.np, fid);
    fread(fp4->ch_z, sizeof(double), g->ab.np, fid);
    fread(fp4->bH_z, sizeof(double), g->ab.np+1, fid);
    fread(fp4->cH_z, sizeof(double), g->ab.np+1, fid);
	
	fp4->k_x  = fp23->k_x  = fp1->k_x;
	fp4->k_z  = fp23->k_z  = fp1->k_z;
	fp4->kh_x = fp23->kh_x = fp1->kh_x;
	fp4->kh_z = fp23->kh_z = fp1->kh_z;
	fp4->kH_x = fp23->kH_x = fp1->kH_x;
	fp4->kH_z = fp23->kH_z = fp1->kH_z;
	
	for ( size_t n=0; n<3; ++n ) {
		fread(mem[n].dx_txx, sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dx_p  , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dx_txz, sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dx_vx , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dx_qx , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dx_vz , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dz_txz, sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fread(mem[n].dz_tzz, sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fread(mem[n].dz_p  , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fread(mem[n].dz_vz , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fread(mem[n].dz_qz , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fread(mem[n].dz_vx , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);		
	}
	
	size_t nnodes = g->nx2 * g->nz2;
    fread(c->mu     , sizeof(double), nnodes, fid);
    fread(c->epsilon, sizeof(double), nnodes, fid);
    fread(c->E      , sizeof(double), nnodes, fid);
    fread(c->M      , sizeof(double), nnodes, fid);
    fread(c->alpha  , sizeof(double), nnodes, fid);
    fread(c->varphi , sizeof(double), nnodes, fid);
    fread(c->tau_s  , sizeof(double), nnodes, fid);
    fread(c->rho_i  , sizeof(double), nnodes, fid);
    fread(c->rho_j  , sizeof(double), nnodes, fid);
    fread(c->rho_f_i, sizeof(double), nnodes, fid);
    fread(c->rho_f_j, sizeof(double), nnodes, fid);
    fread(c->nk_i   , sizeof(double), nnodes, fid);
    fread(c->nk_j   , sizeof(double), nnodes, fid);
    fread(c->mu_ij  , sizeof(double), nnodes, fid);
    fread(c->m_i    , sizeof(double), nnodes, fid);
    fread(c->m_j    , sizeof(double), nnodes, fid);
	
	fread(coeff_v_i, sizeof(double), nnodes, fid);
	fread(coeff_v_j, sizeof(double), nnodes, fid);
	fread(coeff_q_i, sizeof(double), nnodes, fid);
	fread(coeff_q_j, sizeof(double), nnodes, fid);
	fread(W, sizeof(double), 9*nnodes, fid);
	fread(Ws, sizeof(double), 9*nnodes, fid);
	fread(Wtmp, sizeof(double), 9*nnodes, fid);
	fread(delta, sizeof(double), 9*nnodes, fid);
	
	fclose(fid);
}

void save_checkpointVTI(const size_t it, const struct grid *g,
						const struct inputParams *p,
						const struct sourceParams *s,
						const struct outputParams *out,
						const struct saveEnergyVTI *se,
						const struct fac_pml *fp1,
						const struct fac_pml *fp23,
						const struct fac_pml *fp4,
						const struct mem_pml mem[3],
						const struct computationVariablesVTI *c,
						const double *coeff_v_1,
						const double *coeff_v_3,    
						const double *coeff_q_1,
						const double *coeff_q_3,
						const double *W,
						const double *Ws,
						const double *Wtmp,
						const double *delta		 
						) {
	
	struct stat sb;
    if ( stat("chkpt", &sb) < 0 ) {
        if ( errno == ENOENT ) {
            mkdir("chkpt", S_IRWXU);
        }
    }
    char filename[80];
    sprintf(filename, "chkpt/%s_%zd.dat", p->basename, it );
	FILE *fid = fopen(filename, "w");
    if ( fid == NULL ) {
        fprintf(stderr, "Error, cannot open %s for output (%s)\n", filename, strerror(errno));
        exit(1);
    }
	
	fwrite(&it, sizeof(size_t), 1, fid);
	fwrite(g, sizeof(struct grid), 1, fid);
	fwrite(p, sizeof(struct inputParams), 1, fid);
	fwrite(&(s->nsrc), sizeof(size_t), 1, fid);
	fwrite(&(s->nTemplate), sizeof(size_t), 1, fid);
	for ( size_t n=0; n<s->nsrc*s->nTemplate; ++n ) {
		fwrite(&(s->s[n]), sizeof(struct source), 1, fid);
		fwrite(s->s[n].fct, sizeof(double), s->s[n].length, fid);
	}
	fwrite(out, sizeof(struct outputParams), 1, fid);
	for ( size_t nr=0; nr<out->nrec; ++nr ) {
		fwrite(&(out->r[nr]), sizeof(struct record), 1, fid);

        if ( out->r[nr].type == TRACE && p->segy==0 ) {
            size_t nsamples = round( p->duration / out->r[nr].dt );
            fwrite(out->r[nr].data, sizeof(float), nsamples, fid);
		}
        else if ( out->r[nr].type == TRACE && p->segy==0) {
			fflush(out->r[nr].fid);
			rewind(out->r[nr].fid);
			
			if ( stat("chkpt/trc", &sb) < 0 ) {
				if ( errno == ENOENT ) {
					mkdir("chkpt/trc", S_IRWXU);
				}
			}
			char filename2[100];
			sprintf(filename2, "chkpt/trc/%s_%zd", out->basename, it );
			if ( stat(filename2, &sb) < 0 ) {
				if ( errno == ENOENT ) {
					mkdir(filename2, S_IRWXU);
				}
			}		
			sprintf(filename2, "chkpt/trc/%s_%zd/%s_%zd.trc", out->basename, it, out->basename, (nr+1) );
			FILE *fid2 = fopen(filename2, "w");
			if ( fid2 == NULL ) {
				fprintf(stderr, "Error, cannot open %s for output (%s)\n", filename2, strerror(errno));
				exit(1);
			}
			char ch;
			while(!feof(out->r[nr].fid)) {
				ch = fgetc(out->r[nr].fid);
				if(ferror(out->r[nr].fid)) {
					fprintf(stderr, "Error reading source file (trace).\n");
					exit(1);
				}
				if(!feof(out->r[nr].fid)) fputc(ch, fid2);
				if(ferror(fid2)) {
					fprintf(stderr, "Error writing destination file (trace).\n");
					exit(1);
				}
			}
			fclose(fid2);
		}
	}
	if ( p->saveEnergy ) {
		fwrite(se, sizeof(struct saveEnergy), 1, fid);
		size_t npts = (se->i2E-se->i1E+1)*(se->j2E-se->j1E+1);
		fwrite(se->fE1, sizeof(double), npts, fid);
		fwrite(se->fE2, sizeof(double), npts, fid);
		fwrite(se->rEx, sizeof(double), npts, fid);
		fwrite(se->rEz, sizeof(double), npts, fid);
		
		struct stat sb;
		if ( stat("chkpt/E", &sb) < 0 ) {
			if ( errno == ENOENT ) {
				mkdir("chkpt/E", S_IRWXU);
			}
		}
		char fname[100];
		sprintf(fname, "chkpt/E/%s_%zd_E.dat", p->basename, it );
		FILE *fid2 = fopen(fname, "w");
		if ( fid2 == NULL ) {
			fprintf(stderr, "Error, cannot open %s for output (%s)\n", fname, strerror(errno));
			exit(1);
		}
		char ch;
		fflush(se->fid);
		rewind(se->fid);
		while(!feof(se->fid)) {
			ch = fgetc(se->fid);
			if(ferror(se->fid)) {
				fprintf(stderr, "Error reading source file (energy).\n");
				exit(1);
			}
			if(!feof(se->fid)) fputc(ch, fid2);
			if(ferror(fid2)) {
				fprintf(stderr, "Error writing destination file (energy).\n");
				exit(1);
			}
		}
		fclose(fid2);
	}
	
	fwrite(fp1->k_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->k_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->kh_x, sizeof(double), g->ab.np, fid);
	fwrite(fp1->kh_z, sizeof(double), g->ab.np, fid);
    fwrite(fp1->kH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp1->kH_z, sizeof(double), g->ab.np+1, fid);
    fwrite(fp1->b_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->c_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->bh_x, sizeof(double), g->ab.np, fid);
    fwrite(fp1->ch_x, sizeof(double), g->ab.np, fid);
    fwrite(fp1->bH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp1->cH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp1->b_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->c_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp1->bh_z, sizeof(double), g->ab.np, fid);
    fwrite(fp1->ch_z, sizeof(double), g->ab.np, fid);
    fwrite(fp1->bH_z, sizeof(double), g->ab.np+1, fid);
    fwrite(fp1->cH_z, sizeof(double), g->ab.np+1, fid);

    fwrite(fp23->b_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp23->c_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp23->bh_x, sizeof(double), g->ab.np, fid);
    fwrite(fp23->ch_x, sizeof(double), g->ab.np, fid);
    fwrite(fp23->bH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp23->cH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp23->b_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp23->c_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp23->bh_z, sizeof(double), g->ab.np, fid);
    fwrite(fp23->ch_z, sizeof(double), g->ab.np, fid);
    fwrite(fp23->bH_z, sizeof(double), g->ab.np+1, fid);
    fwrite(fp23->cH_z, sizeof(double), g->ab.np+1, fid);

    fwrite(fp4->b_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp4->c_x,  sizeof(double), g->ab.np, fid);
    fwrite(fp4->bh_x, sizeof(double), g->ab.np, fid);
    fwrite(fp4->ch_x, sizeof(double), g->ab.np, fid);
    fwrite(fp4->bH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp4->cH_x, sizeof(double), g->ab.np+1, fid);
    fwrite(fp4->b_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp4->c_z,  sizeof(double), g->ab.np, fid);
    fwrite(fp4->bh_z, sizeof(double), g->ab.np, fid);
    fwrite(fp4->ch_z, sizeof(double), g->ab.np, fid);
    fwrite(fp4->bH_z, sizeof(double), g->ab.np+1, fid);
    fwrite(fp4->cH_z, sizeof(double), g->ab.np+1, fid);
	
	for ( size_t n=0; n<3; ++n ) {
		fwrite(mem[n].dx_txx, sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dx_p  , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dx_txz, sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dx_vx , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dx_qx , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dx_vz , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fwrite(mem[n].dz_txz, sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fwrite(mem[n].dz_tzz, sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fwrite(mem[n].dz_p  , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fwrite(mem[n].dz_vz , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fwrite(mem[n].dz_qz , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fwrite(mem[n].dz_vx , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);		
	}
	
	size_t nnodes = g->nx2 * g->nz2;
	fwrite(c->c11    , sizeof(double), nnodes, fid);
    fwrite(c->c13    , sizeof(double), nnodes, fid);
    fwrite(c->c33    , sizeof(double), nnodes, fid);
    fwrite(c->c55    , sizeof(double), nnodes, fid);
    fwrite(c->epsilon, sizeof(double), nnodes, fid);
    fwrite(c->M      , sizeof(double), nnodes, fid);
    fwrite(c->alpha1 , sizeof(double), nnodes, fid);
    fwrite(c->alpha3 , sizeof(double), nnodes, fid);
    fwrite(c->varphi , sizeof(double), nnodes, fid);
    fwrite(c->tau_s  , sizeof(double), nnodes, fid);
    fwrite(c->rho_i  , sizeof(double), nnodes, fid);
    fwrite(c->rho_j  , sizeof(double), nnodes, fid);
    fwrite(c->rho_f_i, sizeof(double), nnodes, fid);
    fwrite(c->rho_f_j, sizeof(double), nnodes, fid);
    fwrite(c->nk1    , sizeof(double), nnodes, fid);
    fwrite(c->nk3    , sizeof(double), nnodes, fid);
    fwrite(c->m1     , sizeof(double), nnodes, fid);
    fwrite(c->m3     , sizeof(double), nnodes, fid);
	
	fwrite(coeff_v_1, sizeof(double), nnodes, fid);
	fwrite(coeff_v_3, sizeof(double), nnodes, fid);
	fwrite(coeff_q_1, sizeof(double), nnodes, fid);
	fwrite(coeff_q_3, sizeof(double), nnodes, fid);
	fwrite(W, sizeof(double), 9*nnodes, fid);
	fwrite(Ws, sizeof(double), 9*nnodes, fid);
	fwrite(Wtmp, sizeof(double), 9*nnodes, fid);
	fwrite(delta, sizeof(double), 9*nnodes, fid);
	
    if ( p->segy == 1 ) {
        for ( size_t n=0; n<out->nrec; ++n ) {
            if ( out->r[n].type == TRACE ) {
                size_t nsamples = round( p->duration / out->r[n].dt );
                fwrite(out->r[n].data, sizeof(float), nsamples, fid);
            }
        }
    }

	fclose(fid);
}

void read_checkpoint2VTI(const char *filename, const struct grid *g,
						 const struct inputParams *p, struct sourceParams *s,
						 struct outputParams *out, struct saveEnergyVTI *se,
						 struct fac_pml *fp1, struct fac_pml *fp23,
						 struct fac_pml *fp4, struct mem_pml mem[3],
						 struct computationVariablesVTI *c, double *coeff_v_1,
						 double *coeff_v_3, double *coeff_q_1, double *coeff_q_3,
						 double *W, double *Ws, double *Wtmp, double *delta) {
	
	FILE *fid = fopen(filename, "r");
    if ( fid == NULL ) {
        fprintf(stderr, "Error, cannot open %s for input (%s)\n", filename, strerror(errno));
        exit(1);
    }
	size_t it;
	fread(&it, sizeof(size_t), 1, fid);

	fseek(fid, sizeof(struct grid)+sizeof(struct inputParams), SEEK_CUR);
	
	fread(&(s->nsrc), sizeof(size_t), 1, fid);
	fread(&(s->nTemplate), sizeof(size_t), 1, fid);
	if ( NULL == ( s->s = (struct source *) malloc(s->nsrc*s->nTemplate*sizeof(struct source)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	for ( size_t n=0; n<s->nsrc*s->nTemplate; ++n ) {
		fread(&(s->s[n]), sizeof(struct source), 1, fid);
		if ( NULL == ( s->s[n].fct = malloc( s->s[n].length * sizeof(double)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		fread(s->s[n].fct, sizeof(double), s->s[n].length, fid);
	}
	fread(out, sizeof(struct outputParams), 1, fid);
	if ( NULL == ( out->r = (struct record *) malloc(out->nrec*sizeof(struct record)))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
	for ( size_t nr=0; nr<out->nrec; ++nr ) {
		fread(&(out->r[nr]), sizeof(struct record), 1, fid);
		
        if ( out->r[nr].type == TRACE && p->segy == 1 ) {
            size_t nsamples = round( p->duration / out->r[nr].dt );
            if ( NULL == ( out->r[nr].data = (float *) malloc( nsamples * sizeof(float) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
            fread(out->r[nr].data, sizeof(float), nsamples, fid);
        }
        else if ( out->r[nr].type == TRACE && p->segy == 0 ) {
			
            struct stat sb;
            if ( stat("trc", &sb) < 0 ) {
                if ( errno == ENOENT ) {
                    mkdir("trc", S_IRWXU);
                }
            }
            char filename1[100];
            char filename2[80];
            sprintf(filename2, "trc/%s", out->basename );
            if ( stat(filename2, &sb) < 0 ) {
                if ( errno == ENOENT ) {
                    mkdir(filename2, S_IRWXU);
                }
            }
            sprintf(filename1, "chkpt/trc/%s_%zd/%s_%zd.trc", out->basename, it, out->basename, (nr+1) );
            sprintf(filename2, "trc/%s/%s_%zd.trc", out->basename, out->basename, (nr+1) );
            FILE *fid1 = fopen(filename1, "r");
            if ( fid1 == NULL ) {
                fprintf(stderr, "Error, cannot open %s (%s)\n", filename1, strerror(errno));
                exit(1);
            }
            out->r[nr].fid = fopen(filename2, "w+");
            if ( out->r[nr].fid == NULL ) {
                fprintf(stderr, "Error, cannot open %s (%s)\n", filename2, strerror(errno));
                exit(1);
            }
            char ch;
            while(!feof(fid1)) {
                ch = fgetc(fid1);
                if(ferror(fid1)) {
                    fprintf(stderr, "Error reading source file (trace).\n");
                    exit(1);
                }
                if(!feof(fid1)) fputc(ch, out->r[nr].fid);
                if(ferror(out->r[nr].fid)) {
                    fprintf(stderr, "Error writing destination file (trace).\n");
                    exit(1);
                }
            }
            fclose(fid1);
        }
	}
	if ( p->saveEnergy ) {
		fread(se, sizeof(struct saveEnergy), 1, fid);
		size_t npts = (se->i2E-se->i1E+1)*(se->j2E-se->j1E+1);
		if ( NULL == ( se->fE1 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( se->fE2 = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( se->rEx = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		if ( NULL == ( se->rEz = (double *) malloc( npts * sizeof(double) ))) { fprintf(stderr, "Error: cannot allocate memory\n"); abort(); }
		fread(se->fE1, sizeof(double), npts, fid);
		fread(se->fE2, sizeof(double), npts, fid);
		fread(se->rEx, sizeof(double), npts, fid);
		fread(se->rEz, sizeof(double), npts, fid);
		
		struct stat sb;
		if ( stat("E", &sb) < 0 ) {
			if ( errno == ENOENT ) {
				mkdir("E", S_IRWXU);
			}
		}
		char fname1[100];
		char fname[100];
		sprintf(fname1, "chkpt/E/%s_%zd_E.dat", p->basename, it );
		sprintf(fname, "E/%s_E.dat", p->basename );
		FILE *fid1 = fopen(fname1, "r");
		if ( fid1 == NULL ) {
			fprintf(stderr, "Error, cannot open %s (%s)\n", fname1, strerror(errno));
			exit(1);
		}
		se->fid = fopen(fname, "w+");
		if ( se->fid == NULL ) {
			fprintf(stderr, "Error, cannot open %s for output (%s)\n", fname, strerror(errno));
			exit(1);
		}
		char ch;
		while(!feof(fid1)) {
			ch = fgetc(fid1);
			if(ferror(fid1)) {
				fprintf(stderr, "Error reading source file (energy).\n");
				exit(1);
			}
			if(!feof(fid1)) fputc(ch, se->fid);
			if(ferror(se->fid)) {
				fprintf(stderr, "Error writing destination file (energy).\n");
				exit(1);
			}
		}
		fclose(fid1);
	}
	
	fread(fp1->k_x,  sizeof(double), g->ab.np, fid);
    fread(fp1->k_z,  sizeof(double), g->ab.np, fid);
    fread(fp1->kh_x, sizeof(double), g->ab.np, fid);
	fread(fp1->kh_z, sizeof(double), g->ab.np, fid);
    fread(fp1->kH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp1->kH_z, sizeof(double), g->ab.np+1, fid);
    fread(fp1->b_x,  sizeof(double), g->ab.np, fid);
    fread(fp1->c_x,  sizeof(double), g->ab.np, fid);
    fread(fp1->bh_x, sizeof(double), g->ab.np, fid);
    fread(fp1->ch_x, sizeof(double), g->ab.np, fid);
    fread(fp1->bH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp1->cH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp1->b_z,  sizeof(double), g->ab.np, fid);
    fread(fp1->c_z,  sizeof(double), g->ab.np, fid);
    fread(fp1->bh_z, sizeof(double), g->ab.np, fid);
    fread(fp1->ch_z, sizeof(double), g->ab.np, fid);
    fread(fp1->bH_z, sizeof(double), g->ab.np+1, fid);
    fread(fp1->cH_z, sizeof(double), g->ab.np+1, fid);
	
    fread(fp23->b_x,  sizeof(double), g->ab.np, fid);
    fread(fp23->c_x,  sizeof(double), g->ab.np, fid);
    fread(fp23->bh_x, sizeof(double), g->ab.np, fid);
    fread(fp23->ch_x, sizeof(double), g->ab.np, fid);
    fread(fp23->bH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp23->cH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp23->b_z,  sizeof(double), g->ab.np, fid);
    fread(fp23->c_z,  sizeof(double), g->ab.np, fid);
    fread(fp23->bh_z, sizeof(double), g->ab.np, fid);
    fread(fp23->ch_z, sizeof(double), g->ab.np, fid);
    fread(fp23->bH_z, sizeof(double), g->ab.np+1, fid);
    fread(fp23->cH_z, sizeof(double), g->ab.np+1, fid);
	
    fread(fp4->b_x,  sizeof(double), g->ab.np, fid);
    fread(fp4->c_x,  sizeof(double), g->ab.np, fid);
    fread(fp4->bh_x, sizeof(double), g->ab.np, fid);
    fread(fp4->ch_x, sizeof(double), g->ab.np, fid);
    fread(fp4->bH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp4->cH_x, sizeof(double), g->ab.np+1, fid);
    fread(fp4->b_z,  sizeof(double), g->ab.np, fid);
    fread(fp4->c_z,  sizeof(double), g->ab.np, fid);
    fread(fp4->bh_z, sizeof(double), g->ab.np, fid);
    fread(fp4->ch_z, sizeof(double), g->ab.np, fid);
    fread(fp4->bH_z, sizeof(double), g->ab.np+1, fid);
    fread(fp4->cH_z, sizeof(double), g->ab.np+1, fid);

	fp4->k_x  = fp23->k_x  = fp1->k_x;
	fp4->k_z  = fp23->k_z  = fp1->k_z;
	fp4->kh_x = fp23->kh_x = fp1->kh_x;
	fp4->kh_z = fp23->kh_z = fp1->kh_z;
	fp4->kH_x = fp23->kH_x = fp1->kH_x;
	fp4->kH_z = fp23->kH_z = fp1->kH_z;
	
	for ( size_t n=0; n<3; ++n ) {
		fread(mem[n].dx_txx, sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dx_p  , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dx_txz, sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dx_vx , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dx_qx , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dx_vz , sizeof(double), (2*g->ab.np+1)*g->nz2, fid);
		fread(mem[n].dz_txz, sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fread(mem[n].dz_tzz, sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fread(mem[n].dz_p  , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fread(mem[n].dz_vz , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fread(mem[n].dz_qz , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);
		fread(mem[n].dz_vx , sizeof(double), (2*g->ab.np+1)*g->nx2, fid);		
	}
	
	size_t nnodes = g->nx2 * g->nz2;
	fread(c->c11    , sizeof(double), nnodes, fid);
    fread(c->c13    , sizeof(double), nnodes, fid);
    fread(c->c33    , sizeof(double), nnodes, fid);
    fread(c->c55    , sizeof(double), nnodes, fid);
    fread(c->epsilon, sizeof(double), nnodes, fid);
    fread(c->M      , sizeof(double), nnodes, fid);
    fread(c->alpha1 , sizeof(double), nnodes, fid);
    fread(c->alpha3 , sizeof(double), nnodes, fid);
    fread(c->varphi , sizeof(double), nnodes, fid);
    fread(c->tau_s  , sizeof(double), nnodes, fid);
    fread(c->rho_i  , sizeof(double), nnodes, fid);
    fread(c->rho_j  , sizeof(double), nnodes, fid);
    fread(c->rho_f_i, sizeof(double), nnodes, fid);
    fread(c->rho_f_j, sizeof(double), nnodes, fid);
    fread(c->nk1    , sizeof(double), nnodes, fid);
    fread(c->nk3    , sizeof(double), nnodes, fid);
    fread(c->m1     , sizeof(double), nnodes, fid);
    fread(c->m3     , sizeof(double), nnodes, fid);
	
	fread(coeff_v_1, sizeof(double), nnodes, fid);
	fread(coeff_v_3, sizeof(double), nnodes, fid);
	fread(coeff_q_1, sizeof(double), nnodes, fid);
	fread(coeff_q_3, sizeof(double), nnodes, fid);
	fread(W, sizeof(double), 9*nnodes, fid);
	fread(Ws, sizeof(double), 9*nnodes, fid);
	fread(Wtmp, sizeof(double), 9*nnodes, fid);
	fread(delta, sizeof(double), 9*nnodes, fid);
	
	fclose(fid);
}

void save_segy(const struct inputParams *p, const struct outputParams *out,
               const struct sourceParams *src) {
    struct stat sb;
    char fname[80];
    if ( stat("segy", &sb) < 0 ) {
        if ( errno == ENOENT ) {
            mkdir("segy", S_IRWXU);
        }
    }
    sprintf(fname, "segy/%s.sgy", p->basename );
    FILE *fid = fopen(fname, "w");
    if ( fid == NULL ) {
        fprintf(stderr, "Error, cannot open %s for input (%s)\n", fname, strerror(errno));
        exit(1);
    }
    
    if ( verbose > 1 )
        printf("Saving traces in file %s (SEG Y format) ... ", fname);

    /*
     *  Textual file header
     */
    
    char textual_file_header[3200];
    for (size_t n=0; n<3200; ++n) textual_file_header[n] = ' ';
    textual_file_header[0] = 'C';
    textual_file_header[80] = 'C';
//    for (int n=0; n<40; ++n) {
//        sprintf((textual_file_header+n*80), "C% 2d", n+1);
//        textual_file_header[n*80+80] = '\n';
//    }
//    sprintf((textual_file_header+4), "");
//    sprintf((textual_file_header+38*80+4), "SEG Y REV1");
    
    fwrite(&(textual_file_header[0]), sizeof(char), 3200, fid);
    
    
    /*
     *  Binary file header
     */
    
    int32_t i32 = 1;
    fwrite(&i32, 4, 1, fid);  // 1 - job id
    fwrite(&i32, 4, 1, fid);  // 2 - line number
    fwrite(&i32, 4, 1, fid);  // 3 - reel number
    int16_t i16 = 0;
    for (size_t n=0; n<out->nrec; ++n ) {
        if ( out->r[n].type == TRACE ) i16++;
    }
    fwrite(&i16, 2, 1, fid);  // 4 - number of traces
    i16 = 0;
    fwrite(&i16, 2, 1, fid);  // 5 - number of auxiliary traces
    int16_t i16b;
    for (size_t n=0; n<out->nrec; ++n ) {
        if ( out->r[n].type == TRACE ) {
            i16 = 1000000*out->r[n].dt;
            i16b = round( p->duration / out->r[n].dt );
            break;   // we assume all traces have the same dt
        }
    }
    fwrite(&i16, 2, 1, fid);  // 6 - sample int in micros
    fwrite(&i16, 2, 1, fid);  // 7 - orig sample int
    fwrite(&i16b, 2, 1, fid); // 8 - number of samples
    fwrite(&i16b, 2, 1, fid); // 9 - orig number of samples
    i16 = 5;
    fwrite(&i16, 2, 1, fid);  // 10- data in IEEE fp format
    i16 = 0;
    for (size_t n=0; n<out->nrec; ++n ) {
        //  compute fold with Vz component traces
        if ( out->r[n].type == TRACE && out->r[n].comp == VZ ) i16++;
    }
	if ( i16==0 ) {
		// assuming we modelled SH-wave
		for (size_t n=0; n<out->nrec; ++n ) {
			//  compute fold with Vy component traces
			if ( out->r[n].type == TRACE && out->r[n].comp == VY ) i16++;
		}
	}
    fwrite(&i16, 2, 1, fid);  // 11- fold
    i16 = 1;
    fwrite(&i16, 2, 1, fid);  // 12- trace sorting
    fwrite(&i16, 2, 1, fid);  // 13- vertical sum
    i16 = 0;
    fwrite(&i16, 2, 1, fid);  // 14- 
    fwrite(&i16, 2, 1, fid);  // 15- 
    fwrite(&i16, 2, 1, fid);  // 16- 
    fwrite(&i16, 2, 1, fid);  // 17- 
    fwrite(&i16, 2, 1, fid);  // 18- 
    fwrite(&i16, 2, 1, fid);  // 19- 
    fwrite(&i16, 2, 1, fid);  // 20- 
    fwrite(&i16, 2, 1, fid);  // 21- 
    fwrite(&i16, 2, 1, fid);  // 22- 
    fwrite(&i16, 2, 1, fid);  // 23- 
    i16 = 1;
    fwrite(&i16, 2, 1, fid);  // 24- amplitude recovery
    fwrite(&i16, 2, 1, fid);  // 25- meas. system
    fwrite(&i16, 2, 1, fid);  // 26- polarity
    fseek(fid, 3500L, SEEK_SET);
    uint16_t u16 = 1;
    fwrite(&u16, 2, 1, fid);  // 29- SEG Y format revision number
    i16 = 0;
    fwrite(&i16, 2, 1, fid);  // 30- fixed length trace flag
    fwrite(&i16, 2, 1, fid);  // 31- no extended textual file headers
    fseek(fid, 3600L, SEEK_SET);
    
    int32_t tr_no = 1;
    int32_t src_no = p->shotpt_no;
    int32_t ens_no = 1;
    int32_t offset;
    for (size_t n=0; n<out->nrec; ++n ) {
        if ( out->r[n].type == TRACE ) {
            if ( out->r[n].comp != VZ && out->r[n].comp != VY &&
				out->r[n].comp != VX && out->r[n].comp != DIV &&
				out->r[n].comp != CURL ) continue;
            
            /* 
             *  trace header
             */
            fwrite(&tr_no, 4, 1, fid);   // 1 - -> bytes 1-4
            fwrite(&tr_no, 4, 1, fid);   // 2 - -> bytes 5-8
            fwrite(&tr_no, 4, 1, fid);   // 3 - -> bytes 9-12
            fwrite(&tr_no, 4, 1, fid);   // 4 - -> bytes 13-16
            
            fwrite(&src_no, 4, 1, fid);  // 5 - -> bytes 17-20
            if ( p->simulateCMP != 0 ) src_no++;
            fwrite(&ens_no, 4, 1, fid);  // 6 - -> bytes 21-24
            fwrite(&tr_no, 4, 1, fid);   // 7 - -> bytes 25-28
            tr_no++;
            if ( out->r[n].comp == VZ ) i16 = 12;
            else if ( out->r[n].comp == VY ) i16 = 13;
			else if ( out->r[n].comp == VX ) i16 = 14;
            else if ( out->r[n].comp == DIV ) i16 = 23;
            else if ( out->r[n].comp == CURL ) i16 = 24;
            else i16 = 1;
            fwrite(&i16, 2, 1, fid);     // 8 - -> bytes 29-30
            i16 = 1;
            fwrite(&i16, 2, 1, fid);     // 9 - -> bytes 31-32
            fwrite(&i16, 2, 1, fid);     // 10- -> bytes 33-34
            fwrite(&i16, 2, 1, fid);     // 11- -> bytes 35-36
            
            offset = (int32_t)(out->r[n].x - src->s[0].x);
            fwrite(&offset, 4, 1, fid);     // 12- -> bytes 37-40
            i32 = (int32_t)(-out->r[n].z); // rcv grp elevation (z is depth here)
            fwrite(&i32, 4, 1, fid);     // 13- -> bytes 41-44
            i32 = 0;  // (assume surface evel @ source = 0)
            fwrite(&i32, 4, 1, fid);     // 14- -> bytes 45-48
            i32 = (int32_t)(src->s[0].z-i32);
            fwrite(&i32, 4, 1, fid);     // 15- -> bytes 49-52
            i32 = 0;
            fwrite(&i32, 4, 1, fid);     // 16- -> bytes 53-56
            fwrite(&i32, 4, 1, fid);     // 17- -> bytes 57-60
            fwrite(&i32, 4, 1, fid);     // 18- -> bytes 61-64
            fwrite(&i32, 4, 1, fid);     // 19- -> bytes 65-68
            
            i16 = 1;
            fwrite(&i16, 2, 1, fid);     // 20- -> bytes 69-70
            fwrite(&i16, 2, 1, fid);     // 21- -> bytes 71-72
            
            if ( p->simulateCMP == 0 )
                i32 = (int32_t)src->s[0].x;
            else
                i32 = -offset/2;
            fwrite(&i32, 4, 1, fid);     // 22- -> bytes 73-76
            i32 = 0;
            fwrite(&i32, 4, 1, fid);     // 23- -> bytes 77-80
            
            if ( p->simulateCMP == 0 )
                i32 = (int32_t)out->r[n].x;
            else
                i32 = offset/2;
            fwrite(&i32, 4, 1, fid);     // 24- -> bytes 81-84
            i32 = 0;
            fwrite(&i32, 4, 1, fid);     // 25- -> bytes 85-88
            i16 = 1;
            fwrite(&i16, 2, 1, fid);     // 26- -> bytes 89-90
            i16 = 500;   // weathering vel. -> totally arbitrary
            fwrite(&i16, 2, 1, fid);     // 27- -> bytes 91-92
            fwrite(&i16, 2, 1, fid);     // 28- -> bytes 93-94
            i16 = 0;
            fwrite(&i16, 2, 1, fid);     // 29- -> bytes 95-96
            fwrite(&i16, 2, 1, fid);     // 30- -> bytes 97-98
            fwrite(&i16, 2, 1, fid);     // 31- -> bytes 99-100
            fwrite(&i16, 2, 1, fid);     // 32- -> bytes 101-102
            fwrite(&i16, 2, 1, fid);     // 33- -> bytes 103-104
            fwrite(&i16, 2, 1, fid);     // 34- -> bytes 105-106
            fwrite(&i16, 2, 1, fid);     // 35- -> bytes 107-108
            fwrite(&i16, 2, 1, fid);     // 36- -> bytes 109-110
            fwrite(&i16, 2, 1, fid);     // 37- -> bytes 111-112
            fwrite(&i16, 2, 1, fid);     // 38- -> bytes 113-114

            int16_t ns = round( p->duration / out->r[n].dt );
            fwrite(&ns, 2, 1, fid);     // 39- -> bytes 115-116

            i16 = 1000000*out->r[n].dt;
            fwrite(&i16, 2, 1, fid);     // 40- -> bytes 117-118

            i16 = 1;
            fwrite(&i16, 2, 1, fid);     // 41- -> bytes 119-120
            fwrite(&i16, 2, 1, fid);     // 42- -> bytes 121-122
            fwrite(&i16, 2, 1, fid);     // 43- -> bytes 123-124
            
            i16 = 0;
            fwrite(&i16, 2, 1, fid);     // 44- -> bytes 125-126
            fwrite(&i16, 2, 1, fid);     // 45- -> bytes 127-128
            fwrite(&i16, 2, 1, fid);     // 46- -> bytes 129-130
            fwrite(&i16, 2, 1, fid);     // 47- -> bytes 131-132
            fwrite(&i16, 2, 1, fid);     // 48- -> bytes 133-134
            fwrite(&i16, 2, 1, fid);     // 49- -> bytes 135-136
            fwrite(&i16, 2, 1, fid);     // 50- -> bytes 137-138
            fwrite(&i16, 2, 1, fid);     // 51- -> bytes 139-140
            fwrite(&i16, 2, 1, fid);     // 52- -> bytes 141-142
            fwrite(&i16, 2, 1, fid);     // 53- -> bytes 143-144
            fwrite(&i16, 2, 1, fid);     // 54- -> bytes 145-146
            fwrite(&i16, 2, 1, fid);     // 55- -> bytes 147-148
            fwrite(&i16, 2, 1, fid);     // 56- -> bytes 149-150
            fwrite(&i16, 2, 1, fid);     // 57- -> bytes 151-152
            fwrite(&i16, 2, 1, fid);     // 58- -> bytes 153-154
            fwrite(&i16, 2, 1, fid);     // 59- -> bytes 155-156
            
            time_t now;
            struct tm *tm_now;
            
            now = time ( NULL );
            tm_now = gmtime ( &now );
            
            i16 = 1900+tm_now->tm_year;
            fwrite(&i16, 2, 1, fid);     // 60- -> bytes 157-158
            i16 = tm_now->tm_yday;
            fwrite(&i16, 2, 1, fid);     // 61- -> bytes 159-160
            i16 = tm_now->tm_hour;
            fwrite(&i16, 2, 1, fid);     // 62- -> bytes 161-162
            i16 = tm_now->tm_min;
            fwrite(&i16, 2, 1, fid);     // 63- -> bytes 163-164
            i16 = tm_now->tm_sec;
            fwrite(&i16, 2, 1, fid);     // 64- -> bytes 165-166
            i16 = 2;
            fwrite(&i16, 2, 1, fid);     // 65- -> bytes 167-168

            i16 = 0;
            fwrite(&i16, 2, 1, fid);     // 66- -> bytes 169-170
            
            fwrite(&i16, 2, 1, fid);     // 67- -> bytes 171-172
            fwrite(&i16, 2, 1, fid);     // 68- -> bytes 173-174
            fwrite(&i16, 2, 1, fid);     // 69- -> bytes 175-176
            fwrite(&i16, 2, 1, fid);     // 70- -> bytes 177-178

            i16 = 2;
            fwrite(&i16, 2, 1, fid);     // 71- -> bytes 179-180

            i32 = (int32_t)(out->r[n].x - src->s[0].x);
            fwrite(&i32, 4, 1, fid);     // 72- -> bytes 181-184
            i32 = 0;
            fwrite(&i32, 4, 1, fid);     // 73- -> bytes 185-188
            fwrite(&i32, 4, 1, fid);     // 74- -> bytes 189-192
            fwrite(&i32, 4, 1, fid);     // 75- -> bytes 193-196
            fwrite(&i32, 4, 1, fid);     // 76- -> bytes 197-200
            
            i16 = 1;
            fwrite(&i16, 2, 1, fid);     // 77- -> bytes 201-202
            
            i16 = 6;
            fwrite(&i16, 2, 1, fid);     // 78- -> bytes 203-204
            i32 = 1;
            fwrite(&i32, 4, 1, fid);     // 79- -> bytes 205-208
            i16 = 0;
            fwrite(&i16, 2, 1, fid);     // 79- -> bytes 209-210
            i16 = 6;
            fwrite(&i16, 2, 1, fid);     // 80- -> bytes 211-212
            i16 = 0;
            fwrite(&i16, 2, 1, fid);     // 81- -> bytes 213-214
            
            i16 = 1;
            fwrite(&i16, 2, 1, fid);     // 82- -> bytes 215-216
            i16 = 4;
            fwrite(&i16, 2, 1, fid);     // 82- -> bytes 217-218
            
            fseek(fid, 22L, SEEK_CUR);
            
            /*
             *   trace data
             */
            
            fwrite(out->r[n].data, 4, ns, fid);
        }
    }
    
    fclose(fid);
    if ( verbose > 1 )
        printf("done.\n");

}

void check_dt_trc(struct outputParams *out, double dt) {
    
    for (size_t n=0; n<out->nrec; ++n ) {
        if ( out->r[n].type == TRACE ) {
            
            double rem = fmod(1000000.0*out->r[n].dt, 1000000.0*dt );
            if ( rem > 1.e-9 ) {
                fprintf(stdout, "\n *** Warning ***\n");
                fprintf(stdout, "Time step of receiver %zd not a multiple of modeling time step\n",n+1);
                fprintf(stdout, " *** Warning ***\n");
            }
        }
    }
}

