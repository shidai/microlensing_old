// calculate the percentage contribution of optical depth from different types of lens 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

//double l,b;

void display_results (char *title, double result, double error)
{
    printf ("%s ==================\n", title);
    printf ("result = % .6e\n", result);
    printf ("sigma = % .6e\n", error);
    //printf ("error = % .6f = %.2g sigma\n", result - exact,fabs (result - exact) / error);
}

double g (double *k, size_t dim, void *params)
{
    // Ds, l, b: k[0]
    double F,l,b;
    double x0,y0,z0,xs,ys,zs,xb,yb,zb,rs,r0,pmass,p_bulge;
    // Dd: k[1]
    double p_disk,rdd,zdd,p0,beta,ita,h1,h2,H;
    double Deff;
    //double M,mfunc;
    double theta;

    l=k[2];
    b=k[3];

    r0=8.0; // kpc

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc
    //x0=0.89;
    //y0=0.3827;
    //z0=0.2492;  // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    //b=-2.82*3.1415926/180.0;
    //l=1.55*3.1415926/180.0;
    //theta=23.8*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;

    // source distribution: bulge

    pmass=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
    
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    xb=xs*cos(theta)+ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rs=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density

    p_bulge=pmass*exp(-rs);

    // lense distribution: disk
    rdd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
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

    p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);
 
    // mass function of stars
    //if (k[2]>0.7)
    //{
    //    mfunc=pow(k[2]/0.7,-2.0);
    //}
    //else
    //{
    //    mfunc=pow(k[2]/0.7,-1.3);
    //}

    //if (k[2]<=1.0)
    //{
    //    M=k[2];
    //}
    //else if (k[2]>1.0 && k[2]<=8.0)
    //{
    //    M=0.6;
    //}
    //else if (k[2]>8.0 && k[2]<=40.0)
    //{
    //    M=1.35; 
    //}
    //else 
    //{
    //    M=5.0;
    //}    

    if (k[0] > k[1])
    {
        //F=6.029*pow(10.0,-16.0)*k[0]*cos(b)*p_bulge*p_disk*Deff*Deff;
        F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_bulge*p_disk*Deff*Deff;
        //F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_bulge*p_disk*M*mfunc*Deff*Deff/4.3656;
    }
    else
    {
        F=0.0;
    }

    return F;
}

double gb (double *k, size_t dim, void *params)
{
    // Ds, l, b: k[0]
    double F,l,b;
    double x0,y0,z0,xs,ys,zs,xb,yb,zb,rs,r0,pmass,p_bulge;
    // Dd: k[1]
    double xs2,ys2,zs2,xb2,yb2,zb2,rs2,pmass2,p_bulge2;
    double Deff;
    //double M,mfunc;
    double theta;
 
    l=k[2];
    b=k[3];

    r0=8.0; // kpc

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc
    //x0=0.89;
    //y0=0.3827;
    //z0=0.2492;  // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    //b=-2.82*3.1415926/180.0;
    //l=1.55*3.1415926/180.0;
    //theta=23.8*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;

    // source distribution: bulge

    pmass=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
    
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    xb=xs*cos(theta)+ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rs=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density

    p_bulge=pmass*exp(-rs);

    // lenses distribution: bulge

    pmass2=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass2=1.388789*9.6*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
    
    xs2=-k[1]*cos(b)*cos(l)+r0;
    ys2=k[1]*cos(b)*sin(l);
    zs2=k[1]*sin(b);

    xb2=xs2*cos(theta)+ys2*sin(theta);
    yb2=-xs2*sin(theta)+ys2*cos(theta);
    //xb2=xs2*cos(theta)-ys2*sin(theta);
    //yb2=xs2*sin(theta)+ys2*cos(theta);
    zb2=zs2;

    rs2=pow(pow(xb2/x0,2.0)+pow(yb2/y0,2.0)+pow(zb2/z0,2.0),0.5);  // source density

    p_bulge2=pmass2*exp(-rs2);

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);
 
    // mass function of stars
    //if (k[2]>0.7)
    //{
    //    mfunc=pow(k[2]/0.7,-2.0);
    //}
    //else
    //{
    //    mfunc=pow(k[2]/0.7,-1.3);
    //}

    //if (k[2]<=1.0)
    //{
    //    M=k[2];
    //}
    //else if (k[2]>1.0 && k[2]<=8.0)
    //{
    //    M=0.6;
    //}
    //else if (k[2]>8.0 && k[2]<=40.0)
    //{
    //    M=1.35; 
    //}
    //else 
    //{
    //    M=5.0;
    //}    

    if (k[0] > k[1])
    {
        //F=6.029*pow(10.0,-16.0)*k[0]*cos(b)*p_bulge*p_bulge2*Deff*Deff;
        F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_bulge*p_bulge2*Deff*Deff;
        //F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_bulge*p_bulge2*M*mfunc*Deff*Deff/4.3656;
    }
    else
    {
        F=0.0;
    }

    return F;
}

double g_psr_d (double *k, size_t dim, void *params)
{   
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double rd,zd;
    double p_disk,p_pulsar;
    double A,a,B,E;
    double r1;
    double Deff;
    double F;
    double l,b;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;

    l=k[2];
    b=k[3];
    r0=8.0; // kpc

    // source distribution: disk
    rdd=pow(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0,0.5);
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

    p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk

    // lens distribution: psr

    //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
    //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
    //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
    //A=0.58*2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
    //A=1.37*2000.0; // kpc^(-2)    by Kaspi 2006
    A=1.4*2000.0*10000/1.18; // neturon stars, kpc^(-2)    by Kaspi 2006
    //A=2000.0*10000/5.15; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    //B=9.01;
    B=4.01;
    E=0.05; // kpc

    rd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
    zd=k[1]*sin(b);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

    if ( k[0] > k[1] )
    {
        F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_disk*p_pulsar*Deff*Deff;
    }
    else
    {
        F=0.0;
    }

    return F;
}

double source_d (double *k, size_t dim, void *params)
{
    double F,p_disk;
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double l,b;
    //double r0;

    l=k[1];
    b=k[2];
    r0=8.0;
    //l=k[1];
    //b=k[2];

    // source distribution: disk
    rdd=pow(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0,0.5);
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

    p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk

    F=k[0]*k[0]*cos(b)*p_disk;
    //F=k[0]*k[0]*cos(b)*p_bar;

    return F;
}

double g_psr_b (double *k, size_t dim, void *params)
{   
    double rd,zd;
    double x0,y0,z0,rds,pnum;
    double xs,ys,zs,xb,yb,zb;
    double p_pulsar,p_bulge;
    double A,a,B,E;
    double r0,r1;
    double Deff;
    double F;
    double theta;
    double l,b;
 
    l=k[2];
    b=k[3];

    r0=8.0; // kpc

    // source distribution: bulge

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    theta=24.56*3.1415926/180.0;

    pnum=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
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

    rds=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density
    p_bulge=pnum*exp(-rds);

    // lens distribution: psr

    //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
    //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
    //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
    //A=0.58*2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
    //A=1.37*2000.0; // kpc^(-2)    by Kaspi 2006
    A=1.4*2000.0*10000/1.18; // neutron stars, kpc^(-2)    by Kaspi 2006
    //A=2000.0*10000/5.15; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    //B=9.01;
    B=4.01;
    E=0.05; // kpc

    rd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
    zd=k[1]*sin(b);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

    if ( k[0] > k[1] )
    {
        F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_bulge*p_pulsar*Deff*Deff;
    }
    else
    {
        F=0.0;
    }

    return F;
}

double source_b (double *k, size_t dim, void *params)
{
    double F;
    double r0,rds;
    //double r0;
    double x0,y0,z0,pnum,xs,ys,zs;
    //double b,l,xs,ys,zs;
    double theta,xb,yb,zb;
    double l,b;

    l=k[1];
    b=k[2];

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

    rds=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density
    //rds=pow(pow(pow(xb/x0,2.0)+pow(yb/y0,2.0),2.0)+pow(zb/z0,4.0),0.25);  // source density: bulge
    //rds=pow(pow(pow(xs/x0,2.0)+pow(ys/y0,2.0),2.0)+pow(zs/z0,4.0),0.25);  // source density: bulge

    F=k[0]*k[0]*cos(b)*pnum*exp(-rds);
    //F=k[0]*k[0]*cos(b)*p_bar;

    return F;
}

int main (int argc, char *argv[])
{
    double total;
    double bd,ms,wd,bh,ns,ns_p;
    // all stars, disk to bulge
    // positive
    double res,err;

    double xl[4]={3.0,0.0,0.0,-3.1415926/2.0};
    double xu[4]={8.0,8.0,2.0*3.1415926,3.1415926/2.0};
    //double xl[4]={3.0,0.0,0.0,-3.1415926/2.0};
    //double xu[4]={13.0,13.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&g,4,0};

    size_t calls=5000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    // all stars, bulge to bulge
    double res2,err2;

    double xl2[4]={3.0,3.0,0.0,-3.1415926/2.0};
    double xu2[4]={8.0,8.0,2.0*3.1415926,3.1415926/2.0};
    //double xl2[4]={3.0,3.0,0.0,-3.1415926/2.0};
    //double xu2[4]={13.0,13.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T2;
    gsl_rng *r2;

    gsl_monte_function G2={&gb,4,0};

    size_t calls2=5000000;

    T2=gsl_rng_default;

    // neutron stars, pulsar to bulge
    double res_x,err_x;

    double xl_x[4]={3.0,0.0,0.0,-3.1415926/2.0};
    double xu_x[4]={8.0,8.0,2.0*3.1415926,3.1415926/2.0};
    //double xu_x[4]={13.0,13.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_x;
    gsl_rng *r_x;

    gsl_monte_function G_x={&g_psr_b,4,0};

    size_t calls_x=5000000;

    T_x=gsl_rng_default;

    // neutron stars, pulsar to disk
    double res_xd,err_xd;

    double xl_xd[4]={0.0,0.0,0.0,-3.1415926/2.0};
    double xu_xd[4]={8.0,8.0,2.0*3.1415926,3.1415926/2.0};
    //double xu_xd[4]={13.0,13.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_xd;
    gsl_rng *r_xd;

    gsl_monte_function G_xd={&g_psr_d,4,0};

    size_t calls_xd=5000000;

    T_xd=gsl_rng_default;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // source
    double res_s,err_s;

    double xl_s[3]={3.0,0.0,-3.1415926/2.0};
    double xu_s[3]={8.0,2.0*3.1415926,3.1415926/2.0};
    //double xu_s[3]={13.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_s;
    gsl_rng *r_s;

    gsl_monte_function G_s={&source_b,3,0};

    size_t calls_s=5000000;

    T_s=gsl_rng_default;

    r_s = gsl_rng_alloc (T_s);

    gsl_monte_vegas_state *s_s = gsl_monte_vegas_alloc (3);
    gsl_monte_vegas_integrate (&G_s, xl_s, xu_s, 3, 10000, r_s, s_s, &res_s, &err_s);

    do
    {
	gsl_monte_vegas_integrate (&G_s, xl_s, xu_s, 3, calls_s/5, r_s, s_s, &res_s, &err_s);
    }
    while (fabs (gsl_monte_vegas_chisq (s_s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s_s);

    gsl_rng_free (r_s);
    
    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    // source disk  
    double res_sd,err_sd;

    double xl_sd[3]={0.0,0.0,-3.1415926/2.0};
    double xu_sd[3]={8.0,2.0*3.1415926,3.1415926/2.0};
    //double xu_sd[3]={13.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_sd;
    gsl_rng *r_sd;

    gsl_monte_function G_sd={&source_d,3,0};

    size_t calls_sd=50000000;

    T_sd=gsl_rng_default;

    r_sd = gsl_rng_alloc (T_sd);

    gsl_monte_vegas_state *s_sd = gsl_monte_vegas_alloc (3);
    gsl_monte_vegas_integrate (&G_sd, xl_sd, xu_sd, 3, 10000, r_sd, s_sd, &res_sd, &err_sd);

    do
    {
	gsl_monte_vegas_integrate (&G_sd, xl_sd, xu_sd, 3, calls_sd/5, r_sd, s_sd, &res_sd, &err_sd);
    }
    while (fabs (gsl_monte_vegas_chisq (s_sd) - 1.0) > 0.5);

    gsl_monte_vegas_free (s_sd);

    gsl_rng_free (r_sd);
    
    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    // calculation
    
    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (4);
    gsl_monte_vegas_integrate (&G, xl, xu, 4, 10000, r, s,&res, &err);
    //display_results ("vegas warm-up", res, err);
    //printf ("converging...\n");

    do
    {
	gsl_monte_vegas_integrate (&G, xl, xu, 4, calls/5, r, s,&res, &err);
	printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    //display_results ("vegas final", res, err);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r2 = gsl_rng_alloc (T2);

    gsl_monte_vegas_state *s2 = gsl_monte_vegas_alloc (4);
    gsl_monte_vegas_integrate (&G2, xl2, xu2, 4, 10000, r2, s2, &res2, &err2);

    do
    {
	gsl_monte_vegas_integrate (&G2, xl2, xu2, 4, calls2/5, r2, s2, &res2, &err2);
    }
    while (fabs (gsl_monte_vegas_chisq (s2) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r_x=gsl_rng_alloc(T_x);
	
    gsl_monte_vegas_state *s_x = gsl_monte_vegas_alloc (4);
    gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 4, 10000, r_x, s_x,&res_x, &err_x);

    do
    {
	gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 4, calls_x/5, r_x, s_x,&res_x, &err_x);
    }
    while (fabs (gsl_monte_vegas_chisq (s_x) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    r_xd=gsl_rng_alloc(T_xd);
	
    gsl_monte_vegas_state *s_xd = gsl_monte_vegas_alloc (4);
    gsl_monte_vegas_integrate (&G_xd, xl_xd, xu_xd, 4, 10000, r_xd, s_xd, &res_xd, &err_xd);

    do
    {
	gsl_monte_vegas_integrate (&G_xd, xl_xd, xu_xd, 4, calls_xd/5, r_xd, s_xd, &res_xd, &err_xd);
    }
    while (fabs (gsl_monte_vegas_chisq (s_xd) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    total = (fabs(res)/fabs(res_s)+fabs(res2)/fabs(res_s));
    ns_p = (fabs(res_x)+fabs(res_xd))/(fabs(res_s)+fabs(res_sd));
    bd = (0.07*total)/(ns_p+0.94*total);
    ms = (0.62*total)/(ns_p+0.94*total);
    wd = (0.22*total)/(ns_p+0.94*total);
    bh = (0.03*total)/(ns_p+0.94*total);
    ns = ns_p/(ns_p+0.94*total);
    printf ("bd: %e\n", bd);
    printf ("ms: %e\n", ms);
    printf ("wd: %e\n", wd);
    printf ("bh: %e\n", bh);
    printf ("ns: %e\n", ns);
    printf ("ns: %e\n", ns_p);
    printf ("ns: %e\n", 0.06*total);
    printf ("ns: %e %e\n", res_s, res_sd);
    printf ("ns: %e %e\n", 0.06*res/res_s, 0.06*res2/res_s);

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    gsl_monte_vegas_free (s2);

    gsl_rng_free (r2);

    gsl_monte_vegas_free (s_x);

    gsl_rng_free (r_x);

    gsl_monte_vegas_free (s_xd);

    gsl_rng_free (r_xd);

    return 0;
}
