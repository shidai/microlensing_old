// calculate the contribution from disk sources 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

double b=-2.75*3.1415926/180.0;
double l=1.16*3.1415926/180.0;

void display_results (char *title, double result, double error)
{
    printf ("%s ==================\n", title);
    printf ("result = % .6e\n", result);
    printf ("sigma = % .6e\n", error);
    //printf ("error = % .6f = %.2g sigma\n", result - exact,fabs (result - exact) / error);
}

double g1 (double *k, size_t dim, void *params)
{   
    double x0,y0,z0,rds,pnum;
    double xs,ys,zs,xb,yb,zb;
    double p_bulge;
    double r0;
    double Deff;
    double F;
    double theta;

    r0=8.0; // kpc

    // source distribution: bulge

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    theta=24.56*3.1415926/180.0;

    pnum=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    xs=-k[0]*cos(b)*cos(l)+r0;
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    //s=pow(xs*xs+(ys*ys+zs*zs)/(q*q),0.5);
    xb=xs*cos(theta)+ys*sin(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=sqrt(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0));  // source density
    p_bulge=pnum*exp(-rds);

    // lense distribution: bulge
  
    double pnum2,xs2,ys2,zs2,xb2,yb2,zb2,rds2,p_lens;

    pnum2=9.0*pow(10.0,9.0)/1.3888;  // M*kpc^(-3)

    //ys2=k[1]*cos(b)*cos(l)-r0;
    //xs2=-k[1]*cos(b)*sin(l);
    xs2=-k[1]*cos(b)*cos(l)+r0;
    ys2=k[1]*cos(b)*sin(l);
    zs2=k[1]*sin(b);

    xb2=xs2*cos(theta)+ys2*sin(theta);
    //xb2=xs2*cos(theta)-ys2*sin(theta);
    yb2=-xs2*sin(theta)+ys2*cos(theta);
    //yb2=xs2*sin(theta)+ys2*cos(theta);
    zb2=zs2;

    rds2=sqrt(pow(xb2/x0,2.0)+pow(yb2/y0,2.0)+pow(zb2/z0,2.0));  // source density
    p_lens=pnum2*exp(-rds2);

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_bulge*p_lens*Deff*Deff;
    }

    return F;
}

double g2 (double *k, size_t dim, void *params)
{   
    double x0,y0,z0,rds,pnum;
    double xs,ys,zs,xb,yb,zb;
    double p_bulge;
    double r0;
    double Deff;
    double F;
    double theta;

    r0=8.0; // kpc

    // source distribution: bulge

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    theta=24.56*3.1415926/180.0;

    pnum=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    xs=-k[0]*cos(b)*cos(l)+r0;
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    //s=pow(xs*xs+(ys*ys+zs*zs)/(q*q),0.5);
    xb=xs*cos(theta)+ys*sin(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=sqrt(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0));  // source density
    p_bulge=pnum*exp(-rds);

    // lense distribution: disk

    double rdd,zdd,p0,beta,h1,h2,H,ita,p_lens;
    rdd=sqrt(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0);
    zdd=k[1]*sin(b);

    p0=0.0493*pow(10,9); // M*kpc^(-3)
    //p0=1.388789*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=1.3888*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=0.0493*pow(10.0,9.0)/3.1; // kpc^(-3)
    beta=0.565;
    h1=0.270; // kpc
    h2=0.440; // kpc
    H=2.75; // kpc

    if (((rdd/9.025)+0.114)<=0.670)
    {
        ita=0.670;
    }
    else
    {
        ita=(rdd/9.025)+0.114;
    }

    p_lens=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_bulge*p_lens*Deff*Deff;
    }

    return F;
}

double source (double *k, size_t dim, void *params)
{
    double F;
    double r0,rds;
    //double r0;
    double x0,y0,z0,pnum,xs,ys,zs;
    //double b,l,xs,ys,zs;
    double theta,xb,yb,zb;

    // source distribution: bulge
    //x0=1.58;
    //y0=0.62;
    //z0=0.43*8.5/8.0;  // kpc
    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    r0=8.0; // kpc

    //b=0.0;
    //l=0.0;
    //theta=13.4*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;
    //theta=4.56*3.1415926/180.0;

    pnum=9.0*pow(10,9);  // kpc^(-3)
    //pnum=pow(8.5/8.0,3.0)*1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    xb=xs*cos(theta)+ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=sqrt(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0));  // source density
    //rds=pow(pow(pow(xb/x0,2.0)+pow(yb/y0,2.0),2.0)+pow(zb/z0,4.0),0.25);  // source density: bulge
    //rds=pow(pow(pow(xs/x0,2.0)+pow(ys/y0,2.0),2.0)+pow(zs/z0,4.0),0.25);  // source density: bulge

    F=k[0]*k[0]*cos(b)*pnum*exp(-rds);
    //F=k[0]*k[0]*cos(b)*p_bar;

    return F;
}

double g1_d (double *k, size_t dim, void *params)
{   
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double p_disk;
    double Deff;
    double F;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;

    r0=8.0; // kpc

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    // source distribution: disk
    rdd=sqrt(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0);
    zdd=k[0]*sin(b);

    //p0=0.0493*pow(10,9); // M*kpc^(-3)
    p0=1.388789*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=1.3888*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=0.0493*pow(10.0,9.0)/3.1; // kpc^(-3)
    beta=0.565;
    h1=0.270; // kpc
    h2=0.440; // kpc
    H=2.75; // kpc

    if (((rdd/9.025)+0.114)<=0.670)
    {
        ita=0.670;
    }
    else
    {
        ita=(rdd/9.025)+0.114;
    }

    p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)/(cosh(zdd/(ita*h1))*cosh(zdd/(ita*h1)))+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk

    // lense distribution: bulge
  
    double x0,y0,z0,theta,pnum,xs2,ys2,zs2,xb2,yb2,zb2,rds2,p_lens;
    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    theta=24.56*3.1415926/180.0;

    pnum=9.0*pow(10.0,9.0)/1.3888;  // M*kpc^(-3)

    //ys2=k[1]*cos(b)*cos(l)-r0;
    //xs2=-k[1]*cos(b)*sin(l);
    xs2=-k[1]*cos(b)*cos(l)+r0;
    ys2=k[1]*cos(b)*sin(l);
    zs2=k[1]*sin(b);

    xb2=xs2*cos(theta)+ys2*sin(theta);
    //xb2=xs2*cos(theta)-ys2*sin(theta);
    yb2=-xs2*sin(theta)+ys2*cos(theta);
    //yb2=xs2*sin(theta)+ys2*cos(theta);
    zb2=zs2;

    rds2=sqrt(pow(xb2/x0,2.0)+pow(yb2/y0,2.0)+pow(zb2/z0,2.0));  // source density
    p_lens=pnum*exp(-rds2);

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_disk*p_lens*Deff*Deff;
    }

    return F;
}

double g2_d (double *k, size_t dim, void *params)
{   
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double p_disk;
    double Deff;
    double F;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;

    r0=8.0; // kpc

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    // source distribution: disk
    rdd=sqrt(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0);
    zdd=k[0]*sin(b);

    //p0=0.0493*pow(10,9); // M*kpc^(-3)
    p0=1.388789*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=1.3888*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=0.0493*pow(10.0,9.0)/3.1; // kpc^(-3)
    beta=0.565;
    h1=0.270; // kpc
    h2=0.440; // kpc
    H=2.75; // kpc

    if (((rdd/9.025)+0.114)<=0.670)
    {
        ita=0.670;
    }
    else
    {
        ita=(rdd/9.025)+0.114;
    }

    p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)/(cosh(zdd/(ita*h1))*cosh(zdd/(ita*h1)))+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk

    // lense distribution: disk

    double rd,zd,p00,itad,p_lens;
    rd=sqrt(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0);
    zd=k[1]*sin(b);

    p00=0.0493*pow(10,9); // M*kpc^(-3)
    //p0=1.388789*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=1.3888*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=0.0493*pow(10.0,9.0)/3.1; // kpc^(-3)
    beta=0.565;
    h1=0.270; // kpc
    h2=0.440; // kpc
    H=2.75; // kpc

    if (((rd/9.025)+0.114)<=0.670)
    {
        itad=0.670;
    }
    else
    {
        itad=(rd/9.025)+0.114;
    }

    p_lens=(p00/itad)*exp(-(rd-r0)/H)*((1.0-beta)*pow(cosh(zd/(itad*h1)),-2.0)+beta*exp(-fabs(zd)/(itad*h2)));  // :disk

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_disk*p_lens*Deff*Deff;
    }

    return F;
}

double source_d (double *k, size_t dim, void *params)
{
    double F,p_disk;
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    //double r0;

    r0=8.0;
    //l=k[1];
    //b=k[2];

    // source distribution: disk
    rdd=sqrt(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0);
    zdd=k[0]*sin(b);

    //p0=0.0493*pow(10,9); // M*kpc^(-3)
    p0=1.388789*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=1.3888*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=0.0493*pow(10.0,9.0)/3.1; // kpc^(-3)
    beta=0.565;
    h1=0.270; // kpc
    h2=0.440; // kpc
    H=2.75; // kpc

    if (((rdd/9.025)+0.114)<=0.670)
    {
        ita=0.670;
    }
    else
    {
        ita=(rdd/9.025)+0.114;
    }

    p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)/(cosh(zdd/(ita*h1))*cosh(zdd/(ita*h1)))+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk

    F=k[0]*k[0]*cos(b)*p_disk;
    //F=k[0]*k[0]*cos(b)*p_bar;

    return F;
}

int main (int argc, char *argv[])
{
    // source
    double res_s,err_s;

    double xl_s[1]={4.5};
    double xu_s[1]={11.5};

    const gsl_rng_type *T_s;
    gsl_rng *r_s;

    gsl_monte_function G_s={&source,1,0};

    size_t calls_s=5000000;

    gsl_rng_env_setup ();

    T_s=gsl_rng_default;

    // bulge to bulge
    double res1_x,err1_x;

    //double xl_x[4]={3.0,0.0,-1200.0,-1200.0};
    //double xu_x[4]={13.0,13.0,1200.0,1200.0};
    double xl1_x[2]={4.5,4.5};
    double xu1_x[2]={11.5,11.5};

    const gsl_rng_type *T1_x;
    gsl_rng *r1_x;

    gsl_monte_function G1_x={&g1,2,0};

    size_t calls1_x=5000000;

    T1_x=gsl_rng_default;

    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // disk to bulge
    double res2_x,err2_x;

    //double xl_x[4]={3.0,0.0,-1200.0,-1200.0};
    //double xu_x[4]={13.0,13.0,1200.0,1200.0};
    double xl2_x[2]={4.5,0.0};
    double xu2_x[2]={11.5,11.5};

    const gsl_rng_type *T2_x;
    gsl_rng *r2_x;

    gsl_monte_function G2_x={&g2,2,0};

    size_t calls2_x=5000000;

    T2_x=gsl_rng_default;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // source disk
    double res_sd,err_sd;

    double xl_sd[1]={0.0};
    double xu_sd[1]={11.5};

    const gsl_rng_type *T_sd;
    gsl_rng *r_sd;

    gsl_monte_function G_sd={&source_d,1,0};

    size_t calls_sd=5000000;

    T_sd=gsl_rng_default;

    // bulge to disk
    double res1_xd,err1_xd;

    double xl1_xd[2]={0.0,4.5};
    double xu1_xd[2]={11.5,11.5};

    const gsl_rng_type *T1_xd;
    gsl_rng *r1_xd;

    gsl_monte_function G1_xd={&g1_d,2,0};

    size_t calls1_xd=5000000;

    T1_xd=gsl_rng_default;

    ////////////////////////////////////////////////////////////////////////////////////////////

    // disk to disk
    double res2_xd,err2_xd;

    double xl2_xd[2]={0.0,0.0};
    double xu2_xd[2]={11.5,11.5};

    const gsl_rng_type *T2_xd;
    gsl_rng *r2_xd;

    gsl_monte_function G2_xd={&g2_d,2,0};

    size_t calls2_xd=5000000;

    T2_xd=gsl_rng_default;

    ////////////////////////////////////////////////////////////////////////////////////////////

    r_sd = gsl_rng_alloc (T_sd);

    gsl_monte_vegas_state *s_sd = gsl_monte_vegas_alloc (1);

    r1_xd=gsl_rng_alloc(T1_xd);
        
    gsl_monte_vegas_state *s1_xd = gsl_monte_vegas_alloc (2);

    r2_xd=gsl_rng_alloc(T2_xd);
        
    gsl_monte_vegas_state *s2_xd = gsl_monte_vegas_alloc (2);

    r_s = gsl_rng_alloc (T_s);

    gsl_monte_vegas_state *s_s = gsl_monte_vegas_alloc (1);

    r1_x = gsl_rng_alloc(T1_x);
        
    gsl_monte_vegas_state *s1_x = gsl_monte_vegas_alloc (2);

    r2_x = gsl_rng_alloc(T2_x);
        
    gsl_monte_vegas_state *s2_x = gsl_monte_vegas_alloc (2);

    ////////////////////////////////////////////////////////////////////////////////////////

    gsl_monte_vegas_integrate (&G_s, xl_s, xu_s, 1, 10000, r_s, s_s, &res_s, &err_s);

    do
    {
        gsl_monte_vegas_integrate (&G_s, xl_s, xu_s, 1, calls_s/5, r_s, s_s, &res_s, &err_s);
    }
    while (fabs (gsl_monte_vegas_chisq (s_s) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    gsl_monte_vegas_integrate (&G1_x, xl1_x, xu1_x, 2, 10000, r1_x, s1_x,&res1_x, &err1_x);

    do
    {
        gsl_monte_vegas_integrate (&G1_x, xl1_x, xu1_x, 2, calls1_x/5, r1_x, s1_x,&res1_x, &err1_x);
    }
    while (fabs (gsl_monte_vegas_chisq (s1_x) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    gsl_monte_vegas_integrate (&G2_x, xl2_x, xu2_x, 2, 10000, r2_x, s2_x,&res2_x, &err2_x);

    do
    {
        gsl_monte_vegas_integrate (&G2_x, xl2_x, xu2_x, 2, calls2_x/5, r2_x, s2_x,&res2_x, &err2_x);
    }
    while (fabs (gsl_monte_vegas_chisq (s2_x) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    gsl_monte_vegas_integrate (&G_sd, xl_sd, xu_sd, 1, 10000, r_sd, s_sd, &res_sd, &err_sd);

    do
    {
        gsl_monte_vegas_integrate (&G_sd, xl_sd, xu_sd, 1, calls_sd/5, r_sd, s_sd, &res_sd, &err_sd);
    }
    while (fabs (gsl_monte_vegas_chisq (s_sd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    gsl_monte_vegas_integrate (&G1_xd, xl1_xd, xu1_xd, 2, 10000, r1_xd, s1_xd, &res1_xd, &err1_xd);

    do
    {
        gsl_monte_vegas_integrate (&G1_xd, xl1_xd, xu1_xd, 2, calls1_xd/5, r1_xd, s1_xd, &res1_xd, &err1_xd);
    }
    while (fabs (gsl_monte_vegas_chisq (s1_xd) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////


    gsl_monte_vegas_integrate (&G2_xd, xl2_xd, xu2_xd, 2, 10000, r2_xd, s2_xd, &res2_xd, &err2_xd);

    do
    {
        gsl_monte_vegas_integrate (&G2_xd, xl2_xd, xu2_xd, 2, calls2_xd/5, r2_xd, s2_xd, &res2_xd, &err2_xd);
    }
    while (fabs (gsl_monte_vegas_chisq (s2_xd) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    printf ("%e %e %e\n", (fabs(res1_x)+fabs(res2_x))/(fabs(res_s)),(fabs(res1_xd)+fabs(res2_xd))/fabs(res_sd),(fabs(res1_x)+fabs(res2_x)+fabs(res1_xd)+fabs(res2_xd))/(fabs(res_s)+fabs(res_sd)));

    gsl_monte_vegas_free (s1_x);

    gsl_rng_free (r1_x);

    gsl_monte_vegas_free (s2_x);

    gsl_rng_free (r2_x);

    gsl_monte_vegas_free (s_s);

    gsl_rng_free (r_s);

    gsl_monte_vegas_free (s1_xd);

    gsl_rng_free (r1_xd);

    gsl_monte_vegas_free (s2_xd);

    gsl_rng_free (r2_xd);

    gsl_monte_vegas_free (s_sd);

    gsl_rng_free (r_sd);

    fflush (stdout);
    return 0;
}
