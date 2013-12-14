// calculate the percantage contribution of optical depth 
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

double g (double *k, size_t dim, void *params)
{
    double mfunc;

    // mass function of stars
    if (k[0]>0.7)
    {
        mfunc=pow(k[0]/0.7,-2.0);
    }
    else
    {
        mfunc=pow(k[0]/0.7,-1.3);
    }

    return mfunc;
}

double mg (double *k, size_t dim, void *params)
{
    double F,mfunc;

    // mass function of stars
    if (k[0]>0.7)
    {
        mfunc=pow(k[0]/0.7,-2.0);
    }
    else
    {
        mfunc=pow(k[0]/0.7,-1.3);
    }

    F=k[0]*mfunc;

    return F;
}

double g_x (double *k, size_t dim, void *params)
{
    double mfunc,F,M;

    // mass function of stars
    if (k[0]>0.7)
    {
        mfunc=pow(k[0]/0.7,-2.0);
    }
    else
    {
        mfunc=pow(k[0]/0.7,-1.3);
    }

    if (k[0]<=1.0)
    {
        M=k[0];
    }
    else if (k[0]>1.0 && k[0]<=8.0)
    {
	M=0.6;
    }
    else if (k[0]>8.0 && k[0]<=40.0)
    {
    	//M=0.0001;
    	M=1.35;
    }
    else if (k[0]>40.0)
    {
	M=5.0;
    }	

    F=pow(M,1.0)*mfunc;

    return F;
}

int main (void)
{
    //double bd,ms,wd,ns,bh;
    
    // mfunc
    double res,err;

    double xl[1]={0.03};
    double xu[1]={120.0};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&g_x,1,0};

    size_t calls=5000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    // mfunc: ns
    double resn,errn;

    double xln[1]={8.0};
    double xun[1]={40.0};

    const gsl_rng_type *Tn;
    gsl_rng *rn;

    gsl_monte_function Gn={&g,1,0};
    //gsl_monte_function Gn={&g_x,1,0};

    size_t callsn=5000000;

    Tn=gsl_rng_default;

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    // calculation mfunc
    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);
    gsl_monte_vegas_integrate (&G, xl, xu, 1, 10000, r, s, &res, &err);

    do
    {
        gsl_monte_vegas_integrate (&G, xl, xu, 1, calls/5, r, s,&res, &err);
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    rn = gsl_rng_alloc (Tn);

    gsl_monte_vegas_state *sn = gsl_monte_vegas_alloc (1);
    gsl_monte_vegas_integrate (&Gn, xln, xun, 1, 10000, rn, sn, &resn, &errn);

    do
    {
        gsl_monte_vegas_integrate (&Gn, xln, xun, 1, callsn/5, rn, sn, &resn, &errn);
    }
    while (fabs (gsl_monte_vegas_chisq (sn) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    //ns = resn/res;
    printf ("resn: %e\n", resn);
    printf ("res: %e\n", res);
    fflush(stdout);        

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    gsl_monte_vegas_free (sn);

    gsl_rng_free (rn);
    
    return 0;
}
