// normalize potentially observable pulsars 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

void display_results (char *title, double result, double error)
{
    printf ("%s ==================\n", title);
    printf ("result = % .6e\n", result);
    printf ("sigma = % .6e\n", error);
    //printf ("error = % .6f = %.2g sigma\n", result - exact,fabs (result - exact) / error);
}

double g_psr (double *k, size_t dim, void *params)
{   
    double rd,zd;
    double p_pulsar;
    double A,a,B,E;
    double r0,r1;
    double b,l;
    double F;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;

    r0=8.0; // kpc

    l=k[1];
    b=k[2];

    // lens distribution: psr

    //A=1.7*20000.0; // kpc^(-2)    by Kaspi 2006
    A=2.5*1.37*2000.0; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    //B=9.01;
    B=4.01;
    E=0.05; // kpc

    rd=pow(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0,0.5);
    zd=k[0]*sin(b);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density
    //p_pulsar=(A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E))/(k[0]*k[0]);  // pulsar density

    F=k[0]*k[0]*cos(b)*p_pulsar;

    return F;
}

int main (int argc, char *argv[])
{
    double psr;

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    // number of potentially observable pulsars 
    double res_x,err_x;

    double xl_x[3]={0.0,0.0,-3.1415926/2.0};
    double xu_x[3]={8.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_x;
    gsl_rng *r_x;

    gsl_monte_function G_x={&g_psr,3,0};

    size_t calls_x=5000000;

    T_x=gsl_rng_default;

    // calculation
    
    r_x=gsl_rng_alloc(T_x);
        
    gsl_monte_vegas_state *s_x = gsl_monte_vegas_alloc (3);
    gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 3, 10000, r_x, s_x,&res_x, &err_x);

    do
    {
        gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 3, calls_x/5, r_x, s_x,&res_x, &err_x);
    }
    while (fabs (gsl_monte_vegas_chisq (s_x) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    psr = res_x;
    printf ("psr event rate: %e\n", psr);

    gsl_monte_vegas_free (s_x);

    gsl_rng_free (r_x);

    return 0;
}
