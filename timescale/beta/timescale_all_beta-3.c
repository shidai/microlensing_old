// calculate percentage contribution to event rate by different objects
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

double g1 (double *k, size_t dim, void *params)
{
    double F,p_disk;
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double x0,y0,z0,b,l,pnum,xs,ys,zs;
    //double b,l,xs,ys,zs;
    double rds,mfunc;
    //double mfunc;
    //double ifunc,ifunc_vbd;
    double Deff;
    double xl,yl;
    //double v,vl,vb,vrot_disk,vld,vbd,vrot,vbb,vlb;
    // vld,vx,vy,vz are variables: k[3],k[4],k[5],k[6]
    //double vldfunc,vbdfunc,vyfunc,vzfunc;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M,re;
    double delta_func;
    double theta,alpha,xb,yb,zb;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R;

    r0=8.0; // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    b=0.000001;
    l=0.000001;
    //b=5.0*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*sin(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    if ((d*d+R*R-r0*r0)/(2.0*d*R)>1.0)         alpha=acos(1);
    else if ((d*d+R*R-r0*r0)/(2.0*d*R)<-1.0)    alpha=acos(-1);
    else alpha=acos((d*d+R*R-r0*r0)/(2.0*d*R));
    if(sin(l)<0) alpha=-1*alpha;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    // source distribution: bulge
    //x0=1.58;
    //y0=0.62;
    //z0=0.43;  // kpc
    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    //theta=13.4*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;
    //theta=4.56*3.1415926/180.0;

    pnum=9.0*pow(10,9);  // kpc^(-3)
    //pnum=pow(8.5/8.0,3.0)*1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pnum=pow(8.5/8.0,3.0)*0.8549*pow(10.0,9.0);  // M*kpc^(-3)
    //pnum=1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    xs=-k[0]*cos(b)*cos(l)+r0;
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    //s=pow(xs*xs+(ys*ys+zs*zs)/(q*q),0.5);
    xb=xs*cos(theta)+ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=sqrt(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0));  // source density
    //rds=pow(pow(pow(xb/x0,2.0)+pow(yb/y0,2.0),2.0)+pow(zb/z0,4.0),0.25);  // source density: bulge

    //p_bar=pmass*pow(s/abar,-albar)*exp(-(s*s)/(rbar*rbar));

    // lense distribution: disk
    rdd=sqrt(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0);
    zdd=k[1]*sin(b);
    xl=-k[1]*cos(b)*cos(l)+r0;
    yl=k[1]*cos(b)*sin(l);

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

    //if (k[0]<=k[1])
    //{
    //    ifunc=0.0;
    //}
    //else
    //{
    //    ifunc=1.0;
    //}

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    // mass function of stars
    if (k[2]>0.7)
    {
        mfunc=pow(k[2]/0.7,-2.0);
    }
    else
    {
        mfunc=pow(k[2]/0.7,-1.3);
    }

    if (k[2]<=1.0)
    {
        M=k[2];
    }
    else if (k[2]>1.0 && k[2]<=8.0)
    {
        M=0.6;
    }
    else if (k[2]>8.0 && k[2]<=40.0)
    {
        M=0.0; 
    }
    else  
    //else if (k[2]>40.0) 
    {
        M=5.0;
    }    

    //A=0.23;
    //mfunc=fm*k[2]*exp(-A*Deff*Deff*k[2]/(t*t));

    // bulge rotational velocity function
    if ((xs*xs+ys*ys)<1.0)
    {
        vrot=100.0*sqrt(xs*xs+ys*ys); // 100km/s,x/kpc
        //vrot=-100.0*xs; // 100km/s,x/kpc
    }
    else
    {
        vrot=100.0;
        //vrot=-100.0*xs/pow(xs*xs+ys*ys,0.5);
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
    vdx=k[3];
    vdy=k[4];
    vdz=k[5];
    vdl=vrot_disk*cos(alpha)+vdx*sin(l)-vdy*cos(l);
    vdb=-vrot_disk*sin(alpha)*sin(b)+vdx*cos(l)*sin(b)+vdy*sin(l)*sin(b)+vdz*cos(b);
    sigma_vdx=20.0;
    sigma_vdy=30.0;
    sigma_vdz=20.0;

    // source: bulge velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
    vsx=k[6];
    vsy=k[7];
    vsz=k[8];
    vsl=vrot*cos(alpha)+vsx*sin(l)-vsy*cos(l);
    vsb=-vrot*sin(alpha)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
    sigma_vsx=110.0;
    sigma_vsy=82.5;
    sigma_vsz=66.3;

    vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
    vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
    v=sqrt(vl*vl+vb*vb);
    //sigma_y=-5.625*pow(xl*xl+yl*yl,0.5)+75.0;
    //sigma_z=-3.75*pow(xl*xl+yl*yl,0.5)+50.0;
    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    //vlb=-k[4];
    ///////////////////////////////////////////////////////////////////
    //vld=k[3];
    //vbd=k[6];
    //vlb=k[4];
    //vbb=k[5];
    //vl=(vld+vrot_disk-220.0)+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vl=k[3]+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vb=vbd-vbb*(k[1]/k[0]);
    //v=pow(vl*vl+vb*vb,0.5);
    ////////////////////////////////////////////////////////////////////
    //vl=k[3];
    //vb=k[4];
    //mean_vl=(vrot_disk-220.0)+(220.0-vrot)*(k[1]/k[0]);
    //mean_vl=(220.0-vrot)*(k[1]/k[0]);
    //mean_vb=0.0;
    //sigma_vl=pow(sigma_y*sigma_y+82.5*82.5*k[1]*k[1]/(k[0]*k[0]),0.5);
    //sigma_vb=pow(sigma_z*sigma_z+66.3*66.3*k[1]*k[1]/(k[0]*k[0]),0.5);
    //v=pow(vl*vl+vb*vb,0.5);

    //if ((184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl)<=0)
    //{
    //    ifunc_vbd=0.0;
    //}
    //else
    //{
    //    ifunc_vbd=1.0;
    //}

    // vbd as a function of vlb,vx,vy,vz,M
    //vbd=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5)+vbb*(k[1]/k[0]);
    /////////////////////////////////////////////////////////////////////////////////////
    // vb as a function of vl
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    //vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*30.0))*exp(-0.5*vld*vld/900.0);
    //vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*vlb*vlb/(82.5*82.5));
    //vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*vbb*vbb/(66.3*66.3));
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    ////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(vl-mean_vl)*(vl-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));

    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

    re=4.28*pow(10.0,8.0)*Deff*sqrt(M);   // km
    delta_func=(1.0/(sqrt(2.0*3.1415926)*0.01))*exp(-0.5*(log10((re/v)/(24.0*60.0*60.0))-t)*(log10((re/v)/(24.0*60.0*60.0))-t)/(0.01*0.01));

        F=2.77*pow(10.0,-8.0)*(1.0/(k[0]*k[0]*k[0]*k[0]))*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        //F=2.77*pow(10.0,-8.0)*(1.0/(k[0]*k[0]))*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
    }

    return F;
}

double g2 (double *k, size_t dim, void *params)
{
    double F;
    double x0,y0,z0,b,l,pnum,pnum2,xs,ys,zs;
    double xs2,ys2,zs2,r0;
    double rds,rds2,mfunc;
    //double ifunc,ifunc_vbd;
    double Deff;
    //double v,vl,vb,vld,vbd,vlb,vbb,vrot,vrot2;
    //double v,vrot,vrot2;
    // vld,vx,vy,vz are variables: k[3],k[4],k[5],k[6]
    //double vldfunc,vbdfunc,vyfunc,vzfunc;
    //double vl,vb,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M,re;
    double delta_func;
    double theta,alpha,xb,yb,zb,xb2,yb2,zb2;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot2;
    double d,R;

    r0=8.0; // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    b=0.000001;
    l=0.000001;
    //l=5.0*3.1415926/180.0;
    //b=5.0*3.1415926/180.0;

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*sin(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    if ((d*d+R*R-r0*r0)/(2.0*d*R)>1.0)         alpha=acos(1);
    else if ((d*d+R*R-r0*r0)/(2.0*d*R)<-1.0)    alpha=acos(-1);
    else alpha=acos((d*d+R*R-r0*r0)/(2.0*d*R));
    if(sin(l)<0) alpha=-1*alpha;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    // source distribution: bulge
    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc
    //x0=1.58;
    //y0=0.62;
    //z0=0.43;  // kpc

    //theta=13.4*3.1415926/180.0;
    //theta=4.56*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;
    //theta=(90.0-24.56)*3.1415926/180.0;

    pnum=9.0*pow(10.0,9.0);  // kpc^(-3)
    pnum2=1.3888*9.0*pow(10.0,9.0);  // kpc^(-3)
    //pnum2=pow(8.5/8.0,3.0)*0.8549*pow(10.0,9.0);  // M*kpc^(-3)
    //pnum2=pow(8.5/8.0,2.0)*1.3888*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pnum=1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pnum2=1.3888*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pnum2=0.8549*pow(10.0,9.0)/3.1;  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    xb=xs*cos(theta)+ys*sin(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=sqrt(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0));  // source density
    //rds=pow(pow(pow(xb/x0,2.0)+pow(yb/y0,2.0),2.0)+pow(zb/z0,4.0),0.25);  // source density: bulge
    //rds=pow(pow(pow(xs/x0,2.0)+pow(ys/y0,2.0),2.0)+pow(zs/z0,4.0),0.25);  // source density: bulge

    // lense distribution: bulge
  
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
    //rds2=pow(pow(pow(xb2/x0,2.0)+pow(yb2/y0,2.0),2.0)+pow(zb2/z0,4.0),0.25);  // lense density: bulge
    //rds2=pow(pow(pow(xs2/x0,2.0)+pow(ys2/y0,2.0),2.0)+pow(zs2/z0,4.0),0.25);  // lense density: bulge

    //if (k[0]<=k[1])
    //{
    //    ifunc=0.0;
    //}
    //else
    //{
    //    ifunc=1.0;
    //}

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    // mass function of stars
    if (k[2]>0.7)
    {
        mfunc=pow(k[2]/0.7,-2.0);
    }
    else
    {
        mfunc=pow(k[2]/0.7,-1.3);
    }

    if (k[2]<=1.0)
    {
        M=k[2];
    }
    else if (k[2]>1.0 && k[2]<=8.0)
    {
        M=0.6;
        //M=k[2];
    }
    else if (k[2]>8.0 && k[2]<=40.0)
    {
        M=0.0; 
        //M=k[2];
    }
    else 
    {
        M=5.0;
        //M=k[2];
    }    

    //A=0.23;
    //mfunc=fm*k[2]*exp(-A*Deff*Deff*k[2]/(t*t));

    // bulge rotational velocity function: source
    if ((xs*xs+ys*ys)<1.0)
    {
        vrot=100.0*sqrt(xs*xs+ys*ys); // 100km/s,x/kpc
        //vrot=100.0*fabs(xs); // 100km/s,x/kpc
    }
    else
    {
        vrot=100.0;
        //vrot=100.0*fabs(xs)/pow(xs*xs+ys*ys,0.5);
    }

    // bulge rotational velocity function: lens
    if ((xs2*xs2+ys2*ys2)<1.0)
    {
        vrot2=100.0*sqrt(xs2*xs2+ys2*ys2); // 100km/s,x/kpc
        //vrot2=100.0*fabs(xs2); // 100km/s,x/kpc
    }
    else
    {
        vrot2=100.0;
        //vrot2=100.0*fabs(xs2)/pow(xs2*xs2+ys2*ys2,0.5);
    }

    // observor velocity
    vol=-220.0*cos(l);
    vob=220.0*sin(l)*sin(b);

    // deflector: bulge velocity
    // vdx,vdy,vdz,sigma_dx=20.0,sigma_dy=30.0,sigma_dz=20
    vdx=k[3];
    vdy=k[4];
    vdz=k[5];
    vdl=vrot2*cos(alpha)+vdx*sin(l)-vdy*cos(l);
    vdb=-vrot2*sin(alpha)*sin(b)+vdx*cos(l)*sin(b)+vdy*sin(l)*sin(b)+vdz*cos(b);
    sigma_vdx=110.0;
    sigma_vdy=82.5;
    sigma_vdz=66.3;

    // source: bulge velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
    vsx=k[6];
    vsy=k[7];
    vsz=k[8];
    vsl=vrot*cos(alpha)+vsx*sin(l)-vsy*cos(l);
    vsb=-vrot*sin(alpha)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
    sigma_vsx=110.0;
    sigma_vsy=82.5;
    sigma_vsz=66.3;

    vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
    vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
    v=sqrt(vl*vl+vb*vb);
    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    //vlb=-k[4];
    //vlb=k[4];
    //vbb=k[5];
    //vld=k[3];
    //vbd=k[6];
    //vl=(vrot2+vld-220.0)+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vb=vbd-vbb*(k[1]/k[0]);
    //v=pow(vl*vl+vb*vb,0.5);
    
    ///////////////////////////////////////////////////////////////////
    //vl=k[3];
    //vb=k[4];
    //mean_vl=vrot2-vrot*(k[1]/k[0])-220.0*(1.0-k[1]/k[0]);
    //mean_vb=0.0;
    //mean_vb=vrot2-vrot*(k[1]/k[0]);
    //sigma_vl=pow(82.5*82.5*(1.0+k[1]*k[1]/(k[0]*k[0])),0.5);
    //sigma_vb=pow(66.3*66.3*(1.0+k[1]*k[1]/(k[0]*k[0])),0.5);
    //v=pow(vl*vl+vb*vb,0.5);

    //if ((184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl)<=0)
    //{
    //    ifunc_vbd=0.0;
    //}
    //else
    //{
    //    ifunc_vbd=1.0;
    //}

    // vbd as a function of vlb,vx,vy,vz,M
    //vbd=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5)+vbb*(k[1]/k[0]);
    ////////////////////////////////////////////////////////////////////////////////////////////

    // vb as a function of vlb,vx,vy,vz,M
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    //vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*k[3]*k[3]/(82.5*82.5));
    //vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*k[4]*k[4]/(82.5*82.5));
    //vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*k[5]*k[5]/(66.3*66.3));
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*vbd*vbd/(66.3*66.3));
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(vl-mean_vl)*(vl-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));

    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

    re=4.28*pow(10.0,8.0)*Deff*pow(M,0.5);   // km
    // delta_func
    delta_func=(1.0/(sqrt(2.0*3.1415926)*0.01))*exp(-0.5*(log10((re/v)/(24.0*60.0*60.0))-t)*(log10((re/v)/(24.0*60.0*60.0))-t)/(0.01*0.01));

        F=2.77*pow(10.0,-8.0)*(1.0/(k[0]*k[0]*k[0]*k[0]))*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        //F=2.77*pow(10.0,-8.0)*(1.0/(k[0]*k[0]))*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
    }

    return F;
}

double g1d (double *k, size_t dim, void *params)
{
    double F,p_disk;
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double p_disk_s,rd,zd,itad;
    double b,l;
    //double b,l,xs,ys,zs;
    double mfunc;
    //double mfunc;
    //double ifunc,ifunc_vbd;
    double Deff;
    //double v,vl,vb,vrot_disk,vld,vbd,vrot,vbb,vlb;
    // vld,vx,vy,vz are variables: k[3],k[4],k[5],k[6]
    //double vldfunc,vbdfunc,vyfunc,vzfunc;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M,re;
    double delta_func;
    double alpha;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R;

    r0=8.0; // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    b=0.000001;
    l=0.000001;
    //b=5.0*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*sin(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    if ((d*d+R*R-r0*r0)/(2.0*d*R)>1.0)         alpha=acos(1);
    else if ((d*d+R*R-r0*r0)/(2.0*d*R)<-1.0)    alpha=acos(-1);
    else alpha=acos((d*d+R*R-r0*r0)/(2.0*d*R));
    if(sin(l)<0) alpha=-1*alpha;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    // source distribution: disk
    rd=sqrt(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0);
    zd=k[0]*sin(b);

    //p0=0.0493*pow(10,9); // M*kpc^(-3)
    p0=1.388789*0.0493*pow(10.0,9.0); // kpc^(-3)
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

    //p_disk_s=(p0/itad)*exp(-(rd-r0)/H)*((1.0-beta)*pow(cosh(zd/(itad*h1)),-2.0)+beta*exp(-fabs(zd)/(itad*h2)));  // :disk

    if (k[0]<5.0)
    {
        p_disk_s=(p0/itad)*exp(-(rd-r0)/H)*((1.0-beta)*pow(cosh(zd/(itad*h1)),-2.0)+beta*exp(-fabs(zd)/(itad*h2)));  // :disk
    }
    else
    {
        p_disk_s=(1.0/pow(k[0],6.0))*(p0/itad)*exp(-(rd-r0)/H)*((1.0-beta)*pow(cosh(zd/(itad*h1)),-2.0)+beta*exp(-fabs(zd)/(itad*h2)));  // :disk
    }

    // lense distribution: disk
    rdd=sqrt(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0);
    zdd=k[1]*sin(b);

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

    //if (k[0]<=k[1])
    //{
    //    ifunc=0.0;
    //}
    //else
    //{
    //    ifunc=1.0;
    //}

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    // mass function of stars
    if (k[2]>0.7)
    {
        mfunc=pow(k[2]/0.7,-2.0);
    }
    else
    {
        mfunc=pow(k[2]/0.7,-1.3);
    }

    if (k[2]<=1.0)
    {
        M=k[2];
    }
    else if (k[2]>1.0 && k[2]<=8.0)
    {
        M=0.6;
    }
    else if (k[2]>8.0 && k[2]<=40.0)
    {
        M=0.0; 
    }
    else  
    //else if (k[2]>40.0) 
    {
        M=5.0;
    }    

    //A=0.23;
    //mfunc=fm*k[2]*exp(-A*Deff*Deff*k[2]/(t*t));

    // disk rotational velocity function
        vrot=220.0; // 100km/s,x/kpc
        //vrot=-100.0*xs; // 100km/s,x/kpc

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
    vdx=k[3];
    vdy=k[4];
    vdz=k[5];
    vdl=vrot_disk*cos(alpha)+vdx*sin(l)-vdy*cos(l);
    vdb=-vrot_disk*sin(alpha)*sin(b)+vdx*cos(l)*sin(b)+vdy*sin(l)*sin(b)+vdz*cos(b);
    sigma_vdx=20.0;
    sigma_vdy=30.0;
    sigma_vdz=20.0;

    // source: bulge velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
    vsx=k[6];
    vsy=k[7];
    vsz=k[8];
    vsl=vrot*cos(alpha)+vsx*sin(l)-vsy*cos(l);
    vsb=-vrot*sin(alpha)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
    sigma_vsx=20.0;
    sigma_vsy=30.0;
    sigma_vsz=20.0;

    vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
    vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
    v=sqrt(vl*vl+vb*vb);
    //sigma_y=-5.625*pow(xl*xl+yl*yl,0.5)+75.0;
    //sigma_z=-3.75*pow(xl*xl+yl*yl,0.5)+50.0;
    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    //vlb=-k[4];
    ///////////////////////////////////////////////////////////////////
    //vld=k[3];
    //vbd=k[6];
    //vlb=k[4];
    //vbb=k[5];
    //vl=(vld+vrot_disk-220.0)+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vl=k[3]+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vb=vbd-vbb*(k[1]/k[0]);
    //v=pow(vl*vl+vb*vb,0.5);
    ////////////////////////////////////////////////////////////////////
    //vl=k[3];
    //vb=k[4];
    //mean_vl=(vrot_disk-220.0)+(220.0-vrot)*(k[1]/k[0]);
    //mean_vl=(220.0-vrot)*(k[1]/k[0]);
    //mean_vb=0.0;
    //sigma_vl=pow(sigma_y*sigma_y+82.5*82.5*k[1]*k[1]/(k[0]*k[0]),0.5);
    //sigma_vb=pow(sigma_z*sigma_z+66.3*66.3*k[1]*k[1]/(k[0]*k[0]),0.5);
    //v=pow(vl*vl+vb*vb,0.5);

    //if ((184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl)<=0)
    //{
    //    ifunc_vbd=0.0;
    //}
    //else
    //{
    //    ifunc_vbd=1.0;
    //}

    // vbd as a function of vlb,vx,vy,vz,M
    //vbd=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5)+vbb*(k[1]/k[0]);
    /////////////////////////////////////////////////////////////////////////////////////
    // vb as a function of vl
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    //vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*30.0))*exp(-0.5*vld*vld/900.0);
    //vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*vlb*vlb/(82.5*82.5));
    //vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*vbb*vbb/(66.3*66.3));
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    ////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(vl-mean_vl)*(vl-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));

    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

    re=4.28*pow(10.0,8.0)*Deff*sqrt(M);   // km
    delta_func=(1.0/(sqrt(2.0*3.1415926)*0.01))*exp(-0.5*(log10((re/v)/(24.0*60.0*60.0))-t)*(log10((re/v)/(24.0*60.0*60.0))-t)/(0.01*0.01));

        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*p_disk_s*p_disk*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(pow(M,0.5)*mfunc)*(v*vlfunc*vbfunc)*1.02*pow(10.0,-9.0)/4.3656;
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(pow(M,0.5)*mfunc)*(v*vldfunc*vbdfunc*vyfunc*vzfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        //F=2.3*2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(pow(M,-0.5)*mfunc)*(v*vldfunc*vbdfunc*vyfunc*vzfunc)*1.02*pow(10.0,-9.0)/(0.1*4.3656);
        //F=2.3*2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(pow(M,0.5)*mfunc)*(v*vldfunc*vbdfunc*vyfunc*vzfunc)*1.02*pow(10.0,-9.0)/(0.1*4.3656);
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*(Deff*Deff*365.0/t)*(M*mfunc)*(vbfunc*vlfunc);
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-14)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*(Deff*Deff*365.0/t)*(mfunc)*M*(vbfunc*vlfunc)*ifunc*ifunc_vb;
        //F=pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*(Deff*Deff*365.0/t)*(M*mfunc)*(vbdfunc*vldfunc*vyfunc*vzfunc);
    }

    return F;
}

double g2d (double *k, size_t dim, void *params)
{
    double F;
    double p_disk_s,rd,zd;
    double x0,y0,z0,b,l,pnum2,xs2,ys2,zs2;
    double r0;
    double rds2,mfunc;
    //double ifunc,ifunc_vbd;
    double Deff;
    //double v,vl,vb,vld,vbd,vlb,vbb,vrot,vrot2;
    //double v,vrot,vrot2;
    // vld,vx,vy,vz are variables: k[3],k[4],k[5],k[6]
    //double vldfunc,vbdfunc,vyfunc,vzfunc;
    //double vl,vb,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M,re;
    double delta_func;
    double theta,alpha,xb2,yb2,zb2;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot2;
    double d,R;

    r0=8.0; // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    b=0.000001;
    l=0.000001;
    //l=5.0*3.1415926/180.0;
    //b=5.0*3.1415926/180.0;

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*sin(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    if ((d*d+R*R-r0*r0)/(2.0*d*R)>1.0)         alpha=acos(1);
    else if ((d*d+R*R-r0*r0)/(2.0*d*R)<-1.0)    alpha=acos(-1);
    else alpha=acos((d*d+R*R-r0*r0)/(2.0*d*R));
    if(sin(l)<0) alpha=-1*alpha;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    // source distribution: disk
    double p0,beta,h1,h2,H,itad;
    rd=sqrt(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0);
    zd=k[0]*sin(b);

    //p0=0.0493*pow(10,9); // M*kpc^(-3)
    p0=1.388789*0.0493*pow(10.0,9.0); // kpc^(-3)
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

    //p_disk_s=(p0/itad)*exp(-(rd-r0)/H)*((1.0-beta)*pow(cosh(zd/(itad*h1)),-2.0)+beta*exp(-fabs(zd)/(itad*h2)));  // :disk

    if (k[0]<5.0)
    {
        p_disk_s=(p0/itad)*exp(-(rd-r0)/H)*((1.0-beta)*pow(cosh(zd/(itad*h1)),-2.0)+beta*exp(-fabs(zd)/(itad*h2)));  // :disk
    }
    else
    {
        p_disk_s=(1.0/pow(k[0],6.0))*(p0/itad)*exp(-(rd-r0)/H)*((1.0-beta)*pow(cosh(zd/(itad*h1)),-2.0)+beta*exp(-fabs(zd)/(itad*h2)));  // :disk
    }

    // lense distribution: bulge
  
    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc
    //x0=1.58;
    //y0=0.62;
    //z0=0.43;  // kpc

    //theta=13.4*3.1415926/180.0;
    //theta=4.56*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;
    //theta=(90.0-24.56)*3.1415926/180.0;

    pnum2=1.3888*9.0*pow(10.0,9.0);  // kpc^(-3)
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
    //rds2=pow(pow(pow(xb2/x0,2.0)+pow(yb2/y0,2.0),2.0)+pow(zb2/z0,4.0),0.25);  // lense density: bulge
    //rds2=pow(pow(pow(xs2/x0,2.0)+pow(ys2/y0,2.0),2.0)+pow(zs2/z0,4.0),0.25);  // lense density: bulge

    //if (k[0]<=k[1])
    //{
    //    ifunc=0.0;
    //}
    //else
    //{
    //    ifunc=1.0;
    //}

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    // mass function of stars
    if (k[2]>0.7)
    {
        mfunc=pow(k[2]/0.7,-2.0);
    }
    else
    {
        mfunc=pow(k[2]/0.7,-1.3);
    }

    if (k[2]<=1.0)
    {
        M=k[2];
    }
    else if (k[2]>1.0 && k[2]<=8.0)
    {
        M=0.6;
        //M=k[2];
    }
    else if (k[2]>8.0 && k[2]<=40.0)
    {
        M=0.0; 
        //M=k[2];
    }
    else 
    {
        M=5.0;
        //M=k[2];
    }    

    //A=0.23;
    //mfunc=fm*k[2]*exp(-A*Deff*Deff*k[2]/(t*t));

    // disk rotational velocity function: source
        vrot=220.0; // 100km/s,x/kpc
        //vrot=100.0*fabs(xs); // 100km/s,x/kpc

    // bulge rotational velocity function: lens
    if ((xs2*xs2+ys2*ys2)<1.0)
    {
        vrot2=100.0*sqrt(xs2*xs2+ys2*ys2); // 100km/s,x/kpc
        //vrot2=100.0*fabs(xs2); // 100km/s,x/kpc
    }
    else
    {
        vrot2=100.0;
        //vrot2=100.0*fabs(xs2)/pow(xs2*xs2+ys2*ys2,0.5);
    }

    // observor velocity
    vol=-220.0*cos(l);
    vob=220.0*sin(l)*sin(b);

    // deflector: bulge velocity
    // vdx,vdy,vdz,sigma_dx=20.0,sigma_dy=30.0,sigma_dz=20
    vdx=k[3];
    vdy=k[4];
    vdz=k[5];
    vdl=vrot2*cos(alpha)+vdx*sin(l)-vdy*cos(l);
    vdb=-vrot2*sin(alpha)*sin(b)+vdx*cos(l)*sin(b)+vdy*sin(l)*sin(b)+vdz*cos(b);
    sigma_vdx=110.0;
    sigma_vdy=82.5;
    sigma_vdz=66.3;

    // source: bulge velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
    vsx=k[6];
    vsy=k[7];
    vsz=k[8];
    vsl=vrot*cos(alpha)+vsx*sin(l)-vsy*cos(l);
    vsb=-vrot*sin(alpha)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
    sigma_vsx=20.0;
    sigma_vsy=30.0;
    sigma_vsz=20.0;

    vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
    vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
    v=sqrt(vl*vl+vb*vb);
    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    //vlb=-k[4];
    //vlb=k[4];
    //vbb=k[5];
    //vld=k[3];
    //vbd=k[6];
    //vl=(vrot2+vld-220.0)+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vb=vbd-vbb*(k[1]/k[0]);
    //v=pow(vl*vl+vb*vb,0.5);
    
    ///////////////////////////////////////////////////////////////////
    //vl=k[3];
    //vb=k[4];
    //mean_vl=vrot2-vrot*(k[1]/k[0])-220.0*(1.0-k[1]/k[0]);
    //mean_vb=0.0;
    //mean_vb=vrot2-vrot*(k[1]/k[0]);
    //sigma_vl=pow(82.5*82.5*(1.0+k[1]*k[1]/(k[0]*k[0])),0.5);
    //sigma_vb=pow(66.3*66.3*(1.0+k[1]*k[1]/(k[0]*k[0])),0.5);
    //v=pow(vl*vl+vb*vb,0.5);

    //if ((184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl)<=0)
    //{
    //    ifunc_vbd=0.0;
    //}
    //else
    //{
    //    ifunc_vbd=1.0;
    //}

    // vbd as a function of vlb,vx,vy,vz,M
    //vbd=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5)+vbb*(k[1]/k[0]);
    ////////////////////////////////////////////////////////////////////////////////////////////

    // vb as a function of vlb,vx,vy,vz,M
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    //vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*k[3]*k[3]/(82.5*82.5));
    //vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*k[4]*k[4]/(82.5*82.5));
    //vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*k[5]*k[5]/(66.3*66.3));
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*vbd*vbd/(66.3*66.3));
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(vl-mean_vl)*(vl-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));

    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

    re=4.28*pow(10.0,8.0)*Deff*pow(M,0.5);   // km
    // delta_func
    delta_func=(1.0/(sqrt(2.0*3.1415926)*0.01))*exp(-0.5*(log10((re/v)/(24.0*60.0*60.0))-t)*(log10((re/v)/(24.0*60.0*60.0))-t)/(0.01*0.01));

        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*p_disk_s*pnum2*exp(-rds2)*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(pow(M,0.5)*mfunc)*(v*vlfunc*vbfunc)*1.02*pow(10.0,-9.0)/4.3656;
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(pow(M,0.5)*mfunc)*(v*vldfunc*vbdfunc*vyfunc*vzfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        //F=(re/v)*2.3*2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(pow(M,0.5)*mfunc)*(v*vldfunc*vbdfunc*vyfunc*vzfunc)*1.02*pow(10.0,-9.0)/(0.1*4.3656);
        //F=2.3*2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(pow(M,0.5)*mfunc)*(v*vldfunc*vbdfunc*vyfunc*vzfunc)*1.02*pow(10.0,-9.0)/(0.1*4.3656);
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff*Deff*365.0/t)*(M*mfunc)*(vldfunc*vbdfunc*vyfunc*vzfunc)/4.3656;
        //F=2.3*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum2*exp(-0.5*rds2*rds2)*(Deff*Deff*365.0/t)*(mfunc)*M*(vbfunc*vlfunc)*ifunc*ifunc_vb;
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum2*exp(-0.5*rds2*rds2)*(Deff*Deff*365.0/t)*(mfunc)*M*(vbfunc*vlfunc);
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8275*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum2*exp(-0.5*rds2*rds2)*(Deff*Deff*365.0/t)*(mfunc)*(vbdfunc*vldfunc*vyfunc*vzfunc)*ifunc*ifunc_vbd;
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8275*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum2*exp(-0.5*rds2*rds2)*(Deff*Deff*365.0/t)*M*(mfunc)*(vbdfunc*vldfunc*vyfunc*vzfunc);
        //F=pow(2.0*3.1415926,0.5)*110.0*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*1.0726*pow(10,-26)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum*exp(-0.5*rds2*rds2)*ifunc*(Deff*Deff*365.0/t)*(mfunc)*M*(vbdfunc*vldfunc*vyfunc*vzfunc);
        //F=pow(2.0*3.1415926,0.5)*110.0*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*1.0726*pow(10,-26)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum*exp(-0.5*rds2*rds2)*ifunc*ifunc_vbd*(Deff*Deff*365.0/t)*(mfunc)*M*(vbdfunc*vldfunc*vyfunc*vzfunc);
    }

    return F;
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
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M,re;
    double delta_func;
    double theta,alpha;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R;

    r0=8.0; // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    b=0.000001;
    l=0.000001;
    //l=5.0*3.1415926/180.0;
    //b=5.0*3.1415926/180.0;

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*sin(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    if ((d*d+R*R-r0*r0)/(2.0*d*R)>1.0)         alpha=acos(1);
    else if ((d*d+R*R-r0*r0)/(2.0*d*R)<-1.0)    alpha=acos(-1);
    else alpha=acos((d*d+R*R-r0*r0)/(2.0*d*R));
    if(sin(l)<0) alpha=-1*alpha;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    M=1.4;

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

    rds=sqrt(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0));  // source density

    // lens distribution: psr

    //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
    //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
    //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
    //A=0.58*2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
    A=2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
    //A=2000.0*10000/5.15; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    //B=9.01;
    B=4.01;
    E=0.05; // kpc

    rd=sqrt(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0);
    zd=k[1]*sin(b);
    xl=-k[1]*cos(b)*cos(l)+r0;
    yl=k[1]*cos(b)*sin(l);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    // bulge rotational velocity function
    if ((xs*xs+ys*ys)<1.0)
    {
        vrot=100.0*sqrt(xs*xs+ys*ys); // 100km/s,x/kpc
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
    v=sqrt(vl*vl+vb*vb);

    //sigma_y=290.0;
    //sigma_z=290.0;
    //sigma_y=-5.625*pow(xl*xl+yl*yl,0.5)+75.0;
    //sigma_z=-3.75*pow(xl*xl+yl*yl,0.5)+50.0;
    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    //vlb=-k[4];
    ///////////////////////////////////////////////////////////////////
    //vld=k[3];
    //vbd=k[6];
    //vlb=k[4];
    //vbb=k[5];
    //vl=(vld+vrot_disk-220.0)+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vl=k[3]+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vb=vbd-vbb*(k[1]/k[0]);
    //v=pow(vl*vl+vb*vb,0.5);
    ////////////////////////////////////////////////////////////////////
    //vl=k[2];
    //vb=k[3];
    //mean_vl=(vrot_disk-220.0)+(220.0-vrot)*(k[1]/k[0]);
    //mean_vl=(220.0-vrot)*(k[1]/k[0]);
    //mean_vb=0.0;
    //sigma_vl=pow(sigma_y*sigma_y+82.5*82.5*k[1]*k[1]/(k[0]*k[0]),0.5);
    //sigma_vb=pow(sigma_z*sigma_z+66.3*66.3*k[1]*k[1]/(k[0]*k[0]),0.5);
    //v=pow(vl*vl+vb*vb,0.5);

    //if ((184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl)<=0)
    //{
    //    ifunc_vbd=0.0;
    //}
    //else
    //{
    //    ifunc_vbd=1.0;
    //}

    // vbd as a function of vlb,vx,vy,vz,M
    //vbd=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5)+vbb*(k[1]/k[0]);
    /////////////////////////////////////////////////////////////////////////////////////
    // vb as a function of vl
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    //vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*30.0))*exp(-0.5*vld*vld/900.0);
    //vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*vlb*vlb/(82.5*82.5));
    //vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*vbb*vbb/(66.3*66.3));
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    ////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(vl-mean_vl)*(vl-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));

    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

    re=4.28*pow(10.0,8.0)*Deff*sqrt(M);   // km
    // delta_func
    delta_func=(1.0/(sqrt(2.0*3.1415926)*0.01))*exp(-0.5*(log10((re/v)/(24.0*60.0*60.0))-t)*(log10((re/v)/(24.0*60.0*60.0))-t)/(0.01*0.01));

        F=2.77*pow(10.0,-8.0)*(1.0/(k[0]*k[0]*k[0]*k[0]))*cos(b)*pnum*exp(-rds)*p_pulsar*(Deff)*(sqrt(M))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func;
        //F=2.77*pow(10.0,-8.0)*(1.0/(k[0]*k[0]))*cos(b)*pnum*exp(-rds)*p_pulsar*(Deff)*(sqrt(M))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func;
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_pulsar*(Deff)*(sqrt(M))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func;
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

    r0=8.0; // kpc

    //l=k[8];
    //b=k[9];
    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    b=0.000001;
    l=0.000001;
    //b=5.0*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*sin(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    if ((d*d+R*R-r0*r0)/(2.0*d*R)>1.0)         alpha=acos(1);
    else if ((d*d+R*R-r0*r0)/(2.0*d*R)<-1.0)    alpha=acos(-1);
    else alpha=acos((d*d+R*R-r0*r0)/(2.0*d*R));
    if(sin(l)<0) alpha=-1*alpha;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    M=1.4;

    // source distribution: disk
    rdd=sqrt(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0);
    zdd=k[0]*sin(b);
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=k[0]*cos(b)*sin(l);

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

    //p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk

    if (k[0]<5.0)
    {
        p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk
    }
    else
    {
        p_disk=(1.0/pow(k[0],6.0))*(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk
    }

    // lens distribution: psr

    //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
    //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
    //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
    A=2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
    //A=1.37*2000.0; // kpc^(-2)    by Kaspi 2006
    //A=2000.0*10000/5.15; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    //B=9.01;
    B=4.01;
    E=0.05; // kpc

    rd=sqrt(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0);
    zd=k[1]*sin(b);
    xl=-k[1]*cos(b)*cos(l)+r0;
    yl=k[1]*cos(b)*sin(l);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    // disk rotational velocity function
    vrot=220.0;

    // disk rotational velocity function
    //if ((xl*xl+yl*yl)<4.0)
    //{
    //    vrot_disk=220.0*xl; // 100km/s,x/kpc
    //}
    //else
    //{
          vrot_disk=220.0;
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
    v=sqrt(vl*vl+vb*vb);

    //sigma_y=290.0;
    //sigma_z=290.0;
    //sigma_y=-5.625*pow(xl*xl+yl*yl,0.5)+75.0;
    //sigma_z=-3.75*pow(xl*xl+yl*yl,0.5)+50.0;
    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    //vlb=-k[4];
    ///////////////////////////////////////////////////////////////////
    //vld=k[3];
    //vbd=k[6];
    //vlb=k[4];
    //vbb=k[5];
    //vl=(vld+vrot_disk-220.0)+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vl=k[3]+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vb=vbd-vbb*(k[1]/k[0]);
    //v=pow(vl*vl+vb*vb,0.5);
    ////////////////////////////////////////////////////////////////////
    //vl=k[2];
    //vb=k[3];
    //mean_vl=(vrot_disk-220.0)+(220.0-vrot)*(k[1]/k[0]);
    //mean_vl=(220.0-vrot)*(k[1]/k[0]);
    //mean_vb=0.0;
    //sigma_vl=pow(sigma_y*sigma_y+82.5*82.5*k[1]*k[1]/(k[0]*k[0]),0.5);
    //sigma_vb=pow(sigma_z*sigma_z+66.3*66.3*k[1]*k[1]/(k[0]*k[0]),0.5);
    //v=pow(vl*vl+vb*vb,0.5);

    //if ((184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl)<=0)
    //{
    //    ifunc_vbd=0.0;
    //}
    //else
    //{
    //    ifunc_vbd=1.0;
    //}

    // vbd as a function of vlb,vx,vy,vz,M
    //vbd=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5)+vbb*(k[1]/k[0]);
    /////////////////////////////////////////////////////////////////////////////////////
    // vb as a function of vl
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    //vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*30.0))*exp(-0.5*vld*vld/900.0);
    //vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*vlb*vlb/(82.5*82.5));
    //vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*vbb*vbb/(66.3*66.3));
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    ////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(vl-mean_vl)*(vl-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));


    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

    re=4.28*pow(10.0,8.0)*Deff*sqrt(M);   // km
    //ve=365.0*13.57*Deff*pow(M,0.5);

    //vfunc=pow(59.6*Deff*Deff*k[2]/(t*t)-k[3]*k[3],0.5)*exp(0.5*0.002377*k[3]*k[3]);
    // delta_func
    delta_func=(1.0/(sqrt(2.0*3.1415926)*0.01))*exp(-0.5*(log10((re/v)/(24.0*60.0*60.0))-t)*(log10((re/v)/(24.0*60.0*60.0))-t)/(0.01*0.01));

        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*p_disk*p_pulsar*(Deff)*(sqrt(M))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func;
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

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    //b=5.0*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;
    b=0.000001;
    l=0.000001;
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

    F=(1.0/(k[0]*k[0]*k[0]*k[0]))*cos(b)*pnum*exp(-rds);
    //F=(1.0/(k[0]*k[0]))*cos(b)*pnum*exp(-rds);
    //F=k[0]*k[0]*cos(b)*pnum*exp(-rds);
    //F=k[0]*k[0]*cos(b)*p_bar;

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
    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    l=0.000001;
    b=0.000001;
    //b=5.0*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;

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

    //p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk

    if (k[0]<5.0)
    {
        p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk
    }
    else
    {
        p_disk=(1.0/pow(k[0],6.0))*(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk
    }

    F=k[0]*k[0]*cos(b)*p_disk;
    //F=k[0]*k[0]*cos(b)*p_bar;

    return F;
}

int main (int argc, char *argv[])
{
    double galaxy;
    int i;
    int id;  //  process rank
    int p;   //  number of processes

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);
    
    // all stars, disk to bulge
    // positive
    double res,err;

    double xl[9]={3.0,0.0,0.03,-100.0,-150.0,-100.0,-550.0,-450.0,-350.0};
    double xu[9]={13.0,13.0,120.0,100.0,150.0,100.0,550.0,450.0,350.0};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&g1,9,0};

    size_t calls=500000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    // all stars, bulge to bulge
    double res2,err2;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2[9]={3.0,3.0,0.03,-550.0,-450.0,-350.0,-550.0,-450.0,-350.0};
    double xu2[9]={13.0,13.0,120.0,550.0,450.0,350.0,550.0,450.0,350.0};

    const gsl_rng_type *T2;
    gsl_rng *r2;

    gsl_monte_function G2={&g2,9,0};

    size_t calls2=500000000;

    T2=gsl_rng_default;

    // all stars, disk to disk
    // positive
    double resd,errd;

    double xld[9]={0.0,0.0,0.03,-100.0,-150.0,-100.0,-100.0,-150.0,-100.0};
    double xud[9]={13.0,13.0,120.0,100.0,150.0,100.0,100.0,150.0,100.0};

    const gsl_rng_type *Td;
    gsl_rng *rd;

    gsl_monte_function Gd={&g1d,9,0};

    size_t callsd=500000000;

    Td=gsl_rng_default;

    // all stars, bulge to disk
    double res2d,err2d;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2d[9]={0.0,3.0,0.03,-550.0,-450.0,-350.0,-100.0,-150.0,-100.0};
    double xu2d[9]={13.0,13.0,120.0,550.0,450.0,350.0,100.0,150.0,100.0};

    const gsl_rng_type *T2d;
    gsl_rng *r2d;

    gsl_monte_function G2d={&g2d,9,0};

    size_t calls2d=500000000;

    T2d=gsl_rng_default;

    // neutron stars, pulsar to bulge
    double res_x,err_x;

    //double xl_x[4]={3.0,0.0,-1200.0,-1200.0};
    //double xu_x[4]={13.0,13.0,1200.0,1200.0};
    double xl_x[8]={3.0,0.0,-2000.0,-2000.0,-2000.0,-550.0,-450.0,-350.0};
    double xu_x[8]={13.0,13.0,2000.0,2000.0,2000.0,550.0,450.0,350.0};

    const gsl_rng_type *T_x;
    gsl_rng *r_x;

    gsl_monte_function G_x={&g_psr,8,0};

    size_t calls_x=500000000;

    T_x=gsl_rng_default;

    // neutron stars, pulsar to disk 
    double res_xd,err_xd;

    //double xl_x[4]={3.0,0.0,-2000.0,-2000.0};
    //double xu_x[4]={13.0,13.0,2000.0,2000.0};
    double xl_xd[8]={0.0,0.0,-2000.0,-2000.0,-2000.0,-100.0,-150.0,-100.0};
    double xu_xd[8]={13.0,13.0,2000.0,2000.0,2000.0,100.0,150.0,100.0};

    const gsl_rng_type *T_xd;
    gsl_rng *r_xd;

    gsl_monte_function G_xd={&g_psr_disk,8,0};

    size_t calls_xd=500000000;

    T_xd=gsl_rng_default;
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // source
    double res_s,err_s;

    double xl_s[1]={3.0};
    double xu_s[1]={13.0};

    const gsl_rng_type *T_s;
    gsl_rng *r_s;

    gsl_monte_function G_s={&source,1,0};

    size_t calls_s=500000000;

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

    // source disk
    double res_sd,err_sd;

    double xl_sd[1]={0.0};
    double xu_sd[1]={13.0};

    const gsl_rng_type *T_sd;
    gsl_rng *r_sd;

    gsl_monte_function G_sd={&source_disk,1,0};

    size_t calls_sd=500000000;

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
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (9);

    r2 = gsl_rng_alloc (T2);

    gsl_monte_vegas_state *s2 = gsl_monte_vegas_alloc (9);

    rd = gsl_rng_alloc (Td);

    gsl_monte_vegas_state *sd = gsl_monte_vegas_alloc (9);

    r2d = gsl_rng_alloc (T2d);

    gsl_monte_vegas_state *s2d = gsl_monte_vegas_alloc (9);

    r_x=gsl_rng_alloc(T_x);
        
    gsl_monte_vegas_state *s_x = gsl_monte_vegas_alloc (8);

    r_xd=gsl_rng_alloc(T_xd);
        
    gsl_monte_vegas_state *s_xd = gsl_monte_vegas_alloc (8);

    // calculation
    for (i=id;i<=200.0;i+=p)
    {
        //t=pow(10.0,0.015*i);    
        t=0.015*i; 
    
        gsl_monte_vegas_integrate (&G, xl, xu, 9, 10000, r, s,&res, &err);

        do
        {
            gsl_monte_vegas_integrate (&G, xl, xu, 9, calls/5, r, s,&res, &err);
        }
        while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

        gsl_monte_vegas_integrate (&G2, xl2, xu2, 9, 10000, r2, s2, &res2, &err2);

        do
        {
            gsl_monte_vegas_integrate (&G2, xl2, xu2, 9, calls2/5, r2, s2, &res2, &err2);
        }
        while (fabs (gsl_monte_vegas_chisq (s2) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

        gsl_monte_vegas_integrate (&Gd, xld, xud, 9, 10000, rd, sd, &resd, &errd);

        do
        {
            gsl_monte_vegas_integrate (&Gd, xld, xud, 9, callsd/5, rd, sd, &resd, &errd);
        }
        while (fabs (gsl_monte_vegas_chisq (sd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

        gsl_monte_vegas_integrate (&G2d, xl2d, xu2d, 9, 10000, r2d, s2d, &res2d, &err2d);

        do
        {
            gsl_monte_vegas_integrate (&G2d, xl2d, xu2d, 9, calls2d/5, r2d, s2d, &res2d, &err2d);
        }
        while (fabs (gsl_monte_vegas_chisq (s2d) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

        gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 8, 10000, r_x, s_x,&res_x, &err_x);

        do
        {
            gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 8, calls_x/5, r_x, s_x,&res_x, &err_x);
        }
        while (fabs (gsl_monte_vegas_chisq (s_x) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

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

        galaxy = (fabs(res)+fabs(res2)+fabs(resd)+fabs(res2d)+fabs(res_x)+fabs(res_xd))/(fabs(res_s)+fabs(res_sd));
        printf ("%e %e\n", t, galaxy);
        fflush(stdout);        
    }

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    gsl_monte_vegas_free (s2);

    gsl_rng_free (r2);

    gsl_monte_vegas_free (sd);

    gsl_rng_free (rd);

    gsl_monte_vegas_free (s2d);

    gsl_rng_free (r2d);

    gsl_monte_vegas_free (s_x);

    gsl_rng_free (r_x);

    gsl_monte_vegas_free (s_xd);

    gsl_rng_free (r_xd);

    fflush (stdout);
    MPI_Finalize ();
    return 0;
}
