// calculate percentage contribution to event rate by different objects
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
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double rd,zd;
    //double xs,ys;
    //double xs,ys,xl,yl;
    double p_disk,p_pulsar;
    double A,a,B,E;
    double r1;
    double b,l;
    double Deff;
    double F;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M;
    double alpha;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R;

    r0=8.0; // kpc

    l=k[8];
    b=k[9];
    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    //b=0.0;
    //l=0.0;

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {  
    d=k[1]*cos(b);
    //d=k[1]*sin(b);
    R=pow(r0*r0+d*d-2.0*r0*d*cos(l),0.5);
    alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    M=1.4;

    // source distribution: disk
    rdd=pow(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0,0.5);
    zdd=k[0]*sin(b);
    //xs=-k[0]*cos(b)*cos(l)+r0;
    //ys=k[0]*cos(b)*sin(l);

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

    A=1.7*20000.0; // kpc^(-2)    by Kaspi 2006
    //A=1.37*2000.0; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    //B=9.01;
    B=4.01;
    E=0.05; // kpc

    rd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
    zd=k[1]*sin(b);
    //xl=-k[1]*cos(b)*cos(l)+r0;
    //yl=k[1]*cos(b)*sin(l);

    p_pulsar=(A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E))/(k[1]*k[1]);  // pulsar density
    //p_pulsar=(A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E));  // pulsar density

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

    // disk rotational velocity function
    vrot=-220.0;

    // disk rotational velocity function
    //if ((xl*xl+yl*yl)<4.0)
    //{
    //    vrot_disk=220.0*xl; // 100km/s,x/kpc
    //}
    //else
    //{
          vrot_disk=-220.0;
    //}

    // observor velocity
    vol=-220.0*cos(l);
    vob=220.0*sin(l)*sin(b);

    // deflector: disk velocity
    // vdx,vdy,vdz,sigma_dx=20.0,sigma_dy=30.0,sigma_dz=20
    vdx=k[2];
    vdy=k[3];
    vdz=k[4];
    vdl=vrot_disk*cos(alpha)+vdx*sin(l)-vdy*cos(l);
    vdb=-vrot_disk*sin(alpha)*sin(b)+vdx*cos(l)*sin(b)+vdy*sin(l)*sin(b)+vdz*cos(b);
    sigma_vdx=290.0;
    sigma_vdy=290.0;
    sigma_vdz=290.0;

    // source: disk velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
    vsx=k[5];
    vsy=k[6];
    vsz=k[7];
    vsl=vrot*cos(alpha)+vsx*sin(l)-vsy*cos(l);
    vsb=-vrot*sin(alpha)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
    sigma_vsx=20.0;
    sigma_vsy=30.0;
    sigma_vsz=20.0;

    vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
    vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
    v=pow(vl*vl+vb*vb,0.5);


    vdxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*p_disk*p_pulsar*(Deff)*(pow(M,0.5))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0);
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_pulsar*(Deff)*(pow(M,0.5))*(v*vlfunc*vbfunc)*1.02*pow(10.0,-9.0);
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(pow(M,0.5)*mfunc)*(v*vldfunc*vbdfunc*vyfunc*vzfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        //F=2.3*2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(pow(M,-0.5)*mfunc)*(v*vldfunc*vbdfunc*vyfunc*vzfunc)*1.02*pow(10.0,-9.0)/(0.1*4.3656);
        //F=2.3*2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(pow(M,0.5)*mfunc)*(v*vldfunc*vbdfunc*vyfunc*vzfunc)*1.02*pow(10.0,-9.0)/(0.1*4.3656);
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*(Deff*Deff*365.0/t)*(M*mfunc)*(vbfunc*vlfunc);
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-14)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*(Deff*Deff*365.0/t)*(mfunc)*M*(vbfunc*vlfunc)*ifunc*ifunc_vb;
        //F=pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*(Deff*Deff*365.0/t)*(M*mfunc)*(vbdfunc*vldfunc*vyfunc*vzfunc);
    }

    return F;
}

double source (double *k, size_t dim, void *params)
{
    double F,p_disk;
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    //double r0;
    double b,l;

    r0=8.0;
    l=k[1];
    b=k[2];

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

int main (int argc, char *argv[])
{
    double psr;

    // source
    double res_s,err_s;

    //double xl_s[3]={0.0,0.0,0.0};
    //double xu_s[3]={11.5,3.1415926/180.0,3.1415926/180.0};
    double xl_s[3]={0.0,0.0,-3.1415926/2.0};
    double xu_s[3]={11.5,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_s;
    gsl_rng *r_s;

    gsl_monte_function G_s={&source,3,0};

    size_t calls_s=50000000;

    gsl_rng_env_setup ();

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

    // neutron stars, pulsar to disk
    double res_x,err_x;

    //double xl_x[4]={3.0,0.0,-1200.0,-1200.0};
    //double xu_x[4]={13.0,13.0,1200.0,1200.0};
    //double xl_x[10]={0.0,0.01,-2000.0,-2000.0,-2000.0,-100.0,-150.0,-100.0,0.0,0.0};
    //double xu_x[10]={11.5,8.0,2000.0,2000.0,2000.0,100.0,150.0,100.0,3.1415926/180.0,3.1415926/180.0};
    double xl_x[10]={0.0,0.01,-2000.0,-2000.0,-2000.0,-100.0,-150.0,-100.0,0.0,-3.1415926/2.0};
    double xu_x[10]={11.5,8.0,2000.0,2000.0,2000.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_x;
    gsl_rng *r_x;

    gsl_monte_function G_x={&g_psr,10,0};

    size_t calls_x=50000000;

    T_x=gsl_rng_default;

    // calculation
    
    r_x=gsl_rng_alloc(T_x);
        
    gsl_monte_vegas_state *s_x = gsl_monte_vegas_alloc (10);
    gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 10, 10000, r_x, s_x,&res_x, &err_x);

    do
    {
        gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 10, calls_x/5, r_x, s_x,&res_x, &err_x);
    }
    while (fabs (gsl_monte_vegas_chisq (s_x) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    //psr = res_x;
    psr = fabs(res_x)/fabs(res_s);
    //printf ("psr event rate: %e %e %e\n", psr, res_x, res_s);
    printf ("psr event rate: %e\n", psr);

    gsl_monte_vegas_free (s_x);

    gsl_rng_free (r_x);

    return 0;
}
