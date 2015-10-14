//
//  main.c
//  figs_carcione96d
//
//  Created by Bernard Giroux on 2012-10-12.
//
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "V.h"

int main(int argc, const char * argv[])
{
    const double pi = 4.0*atan(1.0);
    
    //
    // sandstone
    //
    double K_s   =   80.0;  // GPa
    double rho_s = 2500.;  // kg/m3
    
    double c11 = 71.8;  // GPa
    double c12 =  3.2;  // GPa
    double c13 =  1.2;  // GPa
    double c33 = 53.4;  // GPa
    double c55 = 26.1;  // GPa
    double phi =  0.2;
    double kappa1 = 600.;  // mD
    double kappa3 = 100.;  // mD
    double T1 = 2.;
    double T3 = 3.6;

    double K_f =  2.5;  // GPa
    double rho_f = 1040.;  // kg/m3
    double eta = 1.;    // cP
    
    
    K_s *= 1000.;                // now in MPa
    K_f *= 1000.;                // now in MPa
    c11  *= 1000.;               // now in MPa
    c12  *= 1000.;               // now in MPa
    c13  *= 1000.;               // now in MPa
    c33  *= 1000.;               // now in MPa
    c55  *= 1000.;               // now in MPa
    eta *= 1.013249965828145e6;  // eta/kappa will be in Mega. kg/m^3/s
    rho_s *= 1.e-6;              // density are now in Gg/m^3 (Mega. kg/m^3)
    rho_f *= 1.e-6;

    
    
    
    
    
    double D = K_s*(1.+phi*(K_s/K_f - 1.));
    
    
    double alpha1 = 1. - (c11+c12+c13)/(3.*K_s);
    double alpha3 = 1. - (2.*c13+c33)/(3.*K_s);
    double M = K_s*K_s/(D-(2.*c11+c33+2.*c12+4.*c13)/9.);
    double rho = (1.-phi)*rho_s + phi*rho_f;
    double m1 = T1*rho_f/phi;
    double m3 = T3*rho_f/phi;
    double omega = 2.*pi*3730.;

    double cu11 = c11 + alpha1*alpha1*M;
    double cu13 = c13 + alpha1*alpha3*M;
    double cu33 = c33 + alpha3*alpha3*M;

    
    double complex * V;
    V = (double complex *)malloc(3*sizeof(double complex));
    
	printf("cu11 = %lg\n", cu11);
	printf("cu13 = %lg\n", cu13);
	printf("cu33 = %lg\n", cu33);
	printf("cu55 = %lg\n", c55);
	printf("alpha1 = %lg\n", alpha1);
	printf("alpha3 = %lg\n", alpha3);
	printf("M = %lg\n", M);
	printf("m1 = %lg\n", m1);
	printf("m3 = %lg\n", m3);
	printf("rho = %lg\n", rho);
	printf("rhof = %lg\n", rho_f);
	printf("omega = %lg\n", omega);
	printf("nk1 = %lg\n", eta/kappa1);
	printf("nk3 = %lg\n", eta/kappa3);

	
    
    FILE * fid = fopen("sandstone.dat","w");
    FILE * fid2 = fopen("sandstone2.dat","w");
    
    for ( double theta=0.0; theta<=pi/2.; theta+= pi/200. ) {
        double lx = cos( theta );
        double lz = sin( theta );
    
        V[0] = VqP1(lx, lz, cu11, cu13, cu33, c55, alpha1, alpha3, M, rho, rho_f, eta/kappa1, eta/kappa3, m1, m3, omega);
        V[1] = VqP2(lx, lz, cu11, cu13, cu33, c55, alpha1, alpha3, M, rho, rho_f, eta/kappa1, eta/kappa3, m1, m3, omega);
        V[2] = VqS(lx, lz, cu11, cu13, cu33, c55, alpha1, alpha3, M, rho, rho_f, eta/kappa1, eta/kappa3, m1, m3, omega);
        
        fprintf(fid, "%lg  %lg  %lg  %lg  %lg  %lg\n", 1000.*lx*creal(1./V[0]),
                1000.*lz*creal(1./V[0]), 1000.*lx*creal(1./V[1]),
                1000.*lz*creal(1./V[1]), 1000.*lx*creal(1./V[2]),
                1000.*lz*creal(1./V[2]));
        
        fprintf(fid2, "%lg  %lg  %lg  %lg  %lg  %lg\n", lx*creal(V[0]),
                lz*creal(V[0]), lx*creal(V[1]),
                lz*creal(V[1]), lx*creal(V[2]),
                lz*creal(V[2]));
        
        

    }
    fclose(fid);
    fclose(fid2);
    
    
    free( V );
    return 0;
}

