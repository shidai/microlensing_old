// calculate event rate and time scale relation of pulsars over bulge star, with pulsar distribution, with monte calo
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <mpi.h>

double t;

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
    double x0,y0,z0,rds,pnum,xs,ys,zs;
    double xb,yb,zb;
    double xl,yl;
    double p_pulsar;
    double A,a,B,E;
    double r0,r1;
    double b,l;
    double Deff;
    double F;
    double M,re;
    double theta;
    double delta_func;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R,alpha;

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
        r0=8.0; // distance from sun to galactic center, kpc

        b=-2.75*3.1415926/180.0;
        l=1.16*3.1415926/180.0;
        //b=0.0;
        //l=0.0;
        //b=5.0*3.1415926/180.0;    // galactic longitude
        //l=5.0*3.1415926/180.0;    // galactic latitude

        // velocity coordinate transformation
        d=k[1]*cos(b);
        //d=k[1]*sin(b);
        R=pow(r0*r0+d*d-2.0*r0*d*cos(l),0.5);
        alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
        //alpha=acos(sqrt(1.0-pow(r0*sin(l)/R,2.0)));
        //alpha=asin(r0*sin(l)/R);
        //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

        // mass of NS
        M=1.4;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // lens distribution: psr

        //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
        //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
        //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
        A=2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
        r1=0.55;  // kpc
        a=1.64;
        //B=8.01;
        B=4.01;
        E=0.05; // kpc

        rd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
        zd=k[1]*sin(b);
        xl=-k[1]*cos(b)*cos(l)+r0;
        yl=k[1]*cos(b)*sin(l);

        p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    
        Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // bulge rotational velocity function
        if ((xs*xs+ys*ys)<1.0)
        {
            vrot=100.0*pow(xs*xs+ys*ys,0.5); // 100km/s,x/kpc
            //vrot=-100.0*xs; // 100km/s,x/kpc
            //vrot=100.0*fabs(xs); // 100km/s,x/kpc
        }
        else
        {
            vrot=100.0;
            //vrot=-100.0*xs/pow(xs*xs+ys*ys,0.5);
            //vrot=100.0*fabs(xs)/pow(xs*xs+ys*ys,0.5);
        }

    // disk rotational velocity function
        //if ((xl*xl+yl*yl)<4.0)
        //{
        //    vrot_disk=220.0*xl; // 100km/s,x/kpc
        //}
        //else
        //{
              vrot_disk=220.0;
              //vrot_disk=-220.0*xl/pow(xl*xl+yl*yl,0.5);
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

    // source: bulge velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
        vsx=k[5];
        vsy=k[6];
        vsz=k[7];
        vsl=vrot*cos(alpha)+vsx*sin(l)-vsy*cos(l);
        vsb=-vrot*sin(alpha)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
        sigma_vsx=110.0;
        sigma_vsy=82.5;
        sigma_vsz=66.3;

        vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
        vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
        v=pow(vl*vl+vb*vb,0.5);

        vdxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
        vdyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
        vdzfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
        vsxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
        vsyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
        vszfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Einstain radius
        re=4.28*pow(10.0,8.0)*Deff*pow(M,0.5);   // km

    // delta_func
        delta_func=(1.0/(pow(2.0*3.1415926,0.5)*0.01))*exp(-0.5*(log10((re/v)/(24.0*60.0*60.0))-t)*(log10((re/v)/(24.0*60.0*60.0))-t)/(0.01*0.01));

        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_pulsar*(Deff)*(pow(M,0.5))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func;
    }

    return F;
}

double source (double *k, size_t dim, void *params)
{
    double F;
    double r0,rds;
    //double r0;
    double x0,y0,z0,b,l,pnum,xs,ys,zs;
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

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;
    //b=5.0*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;
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

double g_psr_disk (double *k, size_t dim, void *params)
{   
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double rd,zd;
    double xs,ys,xl,yl;
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
    double re,delta_func;

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
        r0=8.0; // kpc

        //l=k[8];
        //b=k[9];
        b=-2.75*3.1415926/180.0;
        l=1.16*3.1415926/180.0;
        //b=0.0;
        //l=0.0;
        //b=5.0*3.1415926/180.0;
        //l=5.0*3.1415926/180.0;

        //d=k[1]*sin(b);
        d=k[1]*cos(b);
        R=pow(r0*r0+d*d-2.0*r0*d*cos(l),0.5);
        alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
        //alpha=acos(sqrt(1.0-pow(r0*sin(l)/R,2.0)));
        //alpha=asin(r0*sin(l)/R);

        M=1.4;

    // source distribution: disk
        rdd=pow(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0,0.5);
        zdd=k[0]*sin(b);
        xs=-k[0]*cos(b)*cos(l)+r0;
        ys=k[0]*cos(b)*sin(l);

        //p0=0.0493*pow(10,9); // M*kpc^(-3)
        p0=1.388789*0.0493*pow(10.0,9.0); // kpc^(-3)
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

        A=2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
        r1=0.55;  // kpc
        a=1.64;
        //B=9.01;
        B=4.01;
        E=0.05; // kpc

        rd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
        zd=k[1]*sin(b);
        xl=-k[1]*cos(b)*cos(l)+r0;
        yl=k[1]*cos(b)*sin(l);

        p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

        Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

    // disk rotational velocity function
        vrot=220.0;
        //vrot=-220.0*xs/pow(xs*xs+ys*ys,0.5);

    // disk rotational velocity function
        vrot_disk=220.0;
        //vrot_disk=-220.0*xl/pow(xl*xl+yl*yl,0.5);

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

        re=4.28*pow(10.0,8.0)*Deff*pow(M,0.5);   // km
        //ve=365.0*13.57*Deff*pow(M,0.5);

        //vfunc=pow(59.6*Deff*Deff*k[2]/(t*t)-k[3]*k[3],0.5)*exp(0.5*0.002377*k[3]*k[3]);
        // delta_func
        delta_func=(1.0/(pow(2.0*3.1415926,0.5)*0.01))*exp(-0.5*(log10((re/v)/(24.0*60.0*60.0))-t)*(log10((re/v)/(24.0*60.0*60.0))-t)/(0.01*0.01));

        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*p_disk*p_pulsar*(Deff)*(pow(M,0.5))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func;
    }

    return F;
}

double source_disk (double *k, size_t dim, void *params)
{
    double F,p_disk;
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    //double r0;
    double b,l;

    r0=8.0;
    //l=k[1];
    //b=k[2];
    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;
    //l=0.0;
    //b=0.0;
    //b=5.0*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;

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
    double ratio;
    int i;
    int id;  //  process rank
    int p;   //  number of processes

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);
    
    // neutron stars, pulsar to bulge 
    double res_x,err_x;

    //double xl_x[4]={3.0,0.0,-2000.0,-2000.0};
    //double xu_x[4]={13.0,13.0,2000.0,2000.0};
    double xl_x[8]={3.0,0.0,-20000.0,-20000.0,-20000.0,-550.0,-450.0,-350.0};
    double xu_x[8]={13.0,13.0,20000.0,20000.0,20000.0,550.0,450.0,350.0};

    const gsl_rng_type *T_x;
    gsl_rng *r_x;

    gsl_monte_function G_x={&g_psr,8,0};

    size_t calls_x=5000000;

    gsl_rng_env_setup ();

    T_x=gsl_rng_default;
    
    // source
    double res_s,err_s;

    double xl_s[1]={3.0};
    double xu_s[1]={13.0};

    const gsl_rng_type *T_s;
    gsl_rng *r_s;

    gsl_monte_function G_s={&source,1,0};

    size_t calls_s=5000000;

    T_s=gsl_rng_default;

    r_s = gsl_rng_alloc (T_s);

    gsl_monte_vegas_state *s_s = gsl_monte_vegas_alloc (1);
    gsl_monte_vegas_integrate (&G_s, xl_s, xu_s, 1, 10000, r_s, s_s, &res_s, &err_s);

    do
    {
        gsl_monte_vegas_integrate (&G_s, xl_s, xu_s, 1, calls_s/5, r_s, s_s, &res_s, &err_s);
    }
    while (fabs (gsl_monte_vegas_chisq (s_s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s_s);

    gsl_rng_free (r_s);
    
    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    // neutron stars, pulsar to disk 
    double res_xd,err_xd;

    //double xl_x[4]={3.0,0.0,-2000.0,-2000.0};
    //double xu_x[4]={13.0,13.0,2000.0,2000.0};
    double xl_xd[8]={0.0,0.0,-20000.0,-20000.0,-20000.0,-100.0,-150.0,-100.0};
    double xu_xd[8]={13.0,13.0,20000.0,20000.0,20000.0,100.0,150.0,100.0};

    const gsl_rng_type *T_xd;
    gsl_rng *r_xd;

    gsl_monte_function G_xd={&g_psr_disk,8,0};

    size_t calls_xd=5000000;

    T_xd=gsl_rng_default;
    
    // source disk
    double res_sd,err_sd;

    double xl_sd[1]={0.0};
    double xu_sd[1]={13.0};

    const gsl_rng_type *T_sd;
    gsl_rng *r_sd;

    gsl_monte_function G_sd={&source_disk,1,0};

    size_t calls_sd=5000000;

    T_sd=gsl_rng_default;

    r_sd = gsl_rng_alloc (T_sd);

    gsl_monte_vegas_state *s_sd = gsl_monte_vegas_alloc (1);
    gsl_monte_vegas_integrate (&G_sd, xl_sd, xu_sd, 1, 10000, r_sd, s_sd, &res_sd, &err_sd);

    do
    {
        gsl_monte_vegas_integrate (&G_sd, xl_sd, xu_sd, 1, calls_sd/5, r_sd, s_sd, &res_sd, &err_sd);
    }
    while (fabs (gsl_monte_vegas_chisq (s_sd) - 1.0) > 0.5);

    gsl_monte_vegas_free (s_sd);

    gsl_rng_free (r_sd);
    
    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    // calculation
    for (i=id;i<=200.0;i+=p)
    {
        //t=pow(10.0,0.015*i);    
        t=-0.5+0.0175*i; 

        r_x=gsl_rng_alloc(T_x);
        
        gsl_monte_vegas_state *s_x = gsl_monte_vegas_alloc (8);
        gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 8, 10000, r_x, s_x,&res_x, &err_x);

        // display_results ("vegas warm-up", res, err);
        // printf ("converging...\n");
        do
        {
            gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 8, calls_x/5, r_x, s_x,&res_x, &err_x);
          //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
        }
        while (fabs (gsl_monte_vegas_chisq (s_x) - 1.0) > 0.5);
        // display_results ("vegas final", res, err);

        //////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////

        r_xd=gsl_rng_alloc(T_xd);
        
        gsl_monte_vegas_state *s_xd = gsl_monte_vegas_alloc (8);
        gsl_monte_vegas_integrate (&G_xd, xl_xd, xu_xd, 8, 10000, r_xd, s_xd, &res_xd, &err_xd);

        // display_results ("vegas warm-up", res, err);
        // printf ("converging...\n");
        do
        {
            gsl_monte_vegas_integrate (&G_xd, xl_xd, xu_xd, 8, calls_xd/5, r_xd, s_xd, &res_xd, &err_xd);
          //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
        }
        while (fabs (gsl_monte_vegas_chisq (s_xd) - 1.0) > 0.5);
        // display_results ("vegas final", res, err);

        //////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////

        //ratio = (fabs(res_x)+fabs(res_xd))/(fabs(res_s)+fabs(res_sd));
        ratio = fabs(res_x)/fabs(res_s)+fabs(res_xd)/fabs(res_sd);
        printf ("%e %e\n", t,ratio);
        fflush(stdout);        

        gsl_monte_vegas_free (s_x);

        gsl_rng_free (r_x);
    
        gsl_monte_vegas_free (s_xd);

        gsl_rng_free (r_xd);
    
    }

    //printf ("Process %d is done\n", id);
    fflush (stdout);
    MPI_Finalize ();
    return 0;
}
