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
    //double xl,yl;
    //double v,vl,vb,vrot_disk,vld,vbd,vrot,vbb,vlb;
    // vld,vx,vy,vz are variables: k[3],k[4],k[5],k[6]
    //double vldfunc,vbdfunc,vyfunc,vzfunc;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M;
    double theta,alpha,xb,yb,zb;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R;

    r0=8.0; // kpc

    l=k[9];
    b=k[10];
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
    d=k[1]*sin(b);
    R=pow(r0*r0+d*d-2.0*r0*d*cos(l),0.5);
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

    rds=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density
    //rds=pow(pow(pow(xb/x0,2.0)+pow(yb/y0,2.0),2.0)+pow(zb/z0,4.0),0.25);  // source density: bulge

    //p_bar=pmass*pow(s/abar,-albar)*exp(-(s*s)/(rbar*rbar));

    // lense distribution: disk
    rdd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
    zdd=k[1]*sin(b);
    //xl=-k[1]*cos(b)*cos(l)+r0;
    //yl=k[1]*cos(b)*sin(l);

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

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

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
        vrot=100.0*pow(xs*xs+ys*ys,0.5); // 100km/s,x/kpc
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
    v=pow(vl*vl+vb*vb,0.5);
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

    vdxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

        F=2.77*pow(10.0,-8.0)*(1.0/(k[0]*k[0]))*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(pow(M,0.5)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)/4.3656;
        //F=2.77*pow(10.0,-8.0)*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(pow(M,0.5)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)/4.3656;
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(pow(M,0.5)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)/4.3656;
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
    double M;
    double theta,alpha,xb,yb,zb,xb2,yb2,zb2;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot2;
    double d,R;

    r0=8.0; // kpc

    l=k[9];
    b=k[10];
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
    d=k[1]*sin(b);
    R=pow(r0*r0+d*d-2.0*r0*d*cos(l),0.5);
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

    rds=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density
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

    rds2=pow(pow(xb2/x0,2.0)+pow(yb2/y0,2.0)+pow(zb2/z0,2.0),0.5);  // source density
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

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

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
        vrot=100.0*pow(xs*xs+ys*ys,0.5); // 100km/s,x/kpc
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
        vrot2=100.0*pow(xs2*xs2+ys2*ys2,0.5); // 100km/s,x/kpc
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
    v=pow(vl*vl+vb*vb,0.5);
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

    vdxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

        F=2.77*pow(10.0,-8.0)*(1.0/(k[0]*k[0]))*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(pow(M,0.5)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)/4.3656;
        //F=2.77*pow(10.0,-8.0)*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(pow(M,0.5)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)/4.3656;
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(pow(M,0.5)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)/4.3656;
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
    double M;
    //double delta_func;
    double alpha;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R;

    r0=8.0; // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    //b=0.000001;
    //l=0.000001;
    //b=5.0*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;
    l=k[9];
    b=k[10];

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
        p_disk_s=(p0/itad)*(1.0/(k[0]*k[0]*k[0]*k[0]))*exp(-(rd-r0)/H)*((1.0-beta)*pow(cosh(zd/(itad*h1)),-2.0)+beta*exp(-fabs(zd)/(itad*h2)));  // :disk
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


        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*p_disk_s*p_disk*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)/4.3656;
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
    double M;
    double theta,alpha,xb2,yb2,zb2;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot2;
    double d,R;

    r0=8.0; // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    //b=0.000001;
    //l=0.000001;
    l=k[9];
    b=k[10];

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
        p_disk_s=(p0/itad)*(1.0/(k[0]*k[0]*k[0]*k[0]))*exp(-(rd-r0)/H)*((1.0-beta)*pow(cosh(zd/(itad*h1)),-2.0)+beta*exp(-fabs(zd)/(itad*h2)));  // :disk
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


        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*p_disk_s*pnum2*exp(-rds2)*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)/4.3656;
    }

    return F;
}

double g_psr (double *k, size_t dim, void *params)
{   
    double rd,zd;
    double x0,y0,z0,rds,pnum,xs,ys,zs;
    double xb,yb,zb;
    //double xl,yl;
    double p_pulsar;
    double A,a,B,E;
    double r0,r1;
    double b,l;
    double Deff;
    double F;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M;
    double theta,alpha;
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
    d=k[1]*sin(b);
    R=pow(r0*r0+d*d-2.0*r0*d*cos(l),0.5);
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

    rds=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density

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

    rd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
    zd=k[1]*sin(b);
    //xl=-k[1]*cos(b)*cos(l)+r0;
    //yl=k[1]*cos(b)*sin(l);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

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

    vdxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

        F=2.77*pow(10.0,-8.0)*(1.0/(k[0]*k[0]))*cos(b)*pnum*exp(-rds)*p_pulsar*(Deff)*(pow(M,0.5))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0);
        //F=2.77*pow(10.0,-8.0)*cos(b)*pnum*exp(-rds)*p_pulsar*(Deff)*(pow(M,0.5))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0);
        //F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_pulsar*(Deff)*(pow(M,0.5))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0);
    }

    return F;
}

double g_psr_disk (double *k, size_t dim, void *params)
{   
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double rd,zd;
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
    d=k[1]*sin(b);
    R=pow(r0*r0+d*d-2.0*r0*d*cos(l),0.5);
    if ((d*d+R*R-r0*r0)/(2.0*d*R)>1.0)         alpha=acos(1);
    else if ((d*d+R*R-r0*r0)/(2.0*d*R)<-1.0)    alpha=acos(-1);
    else alpha=acos((d*d+R*R-r0*r0)/(2.0*d*R));
    if(sin(l)<0) alpha=-1*alpha;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
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

    //p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk

    if (k[0]<5.0)
    {
        p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk
    }
    else
    {
        p_disk=(p0/ita)*(1.0/(k[0]*k[0]*k[0]*k[0]))*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk
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

    rd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
    zd=k[1]*sin(b);
    //xl=-k[1]*cos(b)*cos(l)+r0;
    //yl=k[1]*cos(b)*sin(l);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

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
    v=pow(vl*vl+vb*vb,0.5);

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

    vdxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*p_disk*p_pulsar*(Deff)*(pow(M,0.5))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0);
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

    l=k[1];
    b=k[2];
    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
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

    F=(1.0/(k[0]*k[0]))*cos(b)*pnum*exp(-rds);
    //F=cos(b)*pnum*exp(-rds);
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

    //p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk
    if (k[0]<5.0)
    {
        p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk
    }
    else
    {
        p_disk=(p0/ita)*(1.0/(k[0]*k[0]*k[0]*k[0]))*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk
    }

    F=k[0]*k[0]*cos(b)*p_disk;
    //F=k[0]*k[0]*cos(b)*p_bar;

    return F;
}

int main (int argc, char *argv[])
{
    double bd,ms,wd,bh,ns,total;
    double bd_p,ms_p,wd_p,bh_p,ns_p;
    // all stars, disk to bulge
    // positive
    double res,err;

    double xl[11]={3.0,0.0,0.03,-100.0,-150.0,-100.0,-550.0,-450.0,-350,0.0,-3.1415926/2.0};
    double xu[11]={13.0,13.0,120.0,100.0,150.0,100.0,550.0,450.0,350.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&g1,11,0};

    size_t calls=50000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    // all stars, bulge to bulge
    double res2,err2;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2[11]={3.0,3.0,0.03,-550.0,-450.0,-350.0,-550.0,-450.0,-350,0.0,-3.1415926/2.0};
    double xu2[11]={13.0,13.0,120.0,550.0,450.0,350.0,550.0,450.0,350.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T2;
    gsl_rng *r2;

    gsl_monte_function G2={&g2,11,0};

    size_t calls2=50000000;

    T2=gsl_rng_default;

    // all stars, disk to disk
    // positive
    double resd,errd;

    double xld[11]={0.0,0.0,0.03,-100.0,-150.0,-100.0,-100.0,-150.0,-100,0.0,-3.1415926/2.0};
    double xud[11]={13.0,13.0,120.0,100.0,150.0,100.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *Td;
    gsl_rng *rd;

    gsl_monte_function Gd={&g1d,11,0};

    size_t callsd=50000000;

    Td=gsl_rng_default;

    // all stars, bulge to disk 
    double res2d,err2d;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2d[11]={0.0,3.0,0.03,-550.0,-450.0,-350.0,-100.0,-150.0,-100,0.0,-3.1415926/2.0};
    double xu2d[11]={13.0,13.0,120.0,550.0,450.0,350.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T2d;
    gsl_rng *r2d;

    gsl_monte_function G2d={&g2d,11,0};

    size_t calls2d=50000000;

    T2d=gsl_rng_default;

    // neutron stars, pulsar to bulge
    double res_x,err_x;

    //double xl_x[4]={3.0,0.0,-1200.0,-1200.0};
    //double xu_x[4]={13.0,13.0,1200.0,1200.0};
    double xl_x[10]={3.0,0.0,-2000.0,-2000.0,-2000.0,-550.0,-450.0,-350,0.0,-3.1415926/2.0};
    double xu_x[10]={13.0,13.0,2000.0,2000.0,2000.0,550.0,450.0,350.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_x;
    gsl_rng *r_x;

    gsl_monte_function G_x={&g_psr,10,0};

    size_t calls_x=50000000;

    T_x=gsl_rng_default;

    // neutron stars, pulsar to disk
    double res_xd,err_xd;

    //double xl_x[4]={3.0,0.0,-1200.0,-1200.0};
    //double xu_x[4]={13.0,13.0,1200.0,1200.0};
    double xl_xd[10]={0.0,0.0,-2000.0,-2000.0,-2000.0,-100.0,-150.0,-100,0.0,-3.1415926/2.0};
    double xu_xd[10]={13.0,13.0,2000.0,2000.0,2000.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_xd;
    gsl_rng *r_xd;

    gsl_monte_function G_xd={&g_psr_disk,10,0};

    size_t calls_xd=50000000;

    T_xd=gsl_rng_default;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // bd, disk to bulge
    double res_bd,err_bd;

    double xl_bd[11]={3.0,0.0,0.03,-100.0,-150.0,-100.0,-550.0,-450.0,-350,0.0,-3.1415926/2.0};
    double xu_bd[11]={13.0,13.0,0.08,100.0,150.0,100.0,550.0,450.0,350.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_bd;
    gsl_rng *r_bd;

    gsl_monte_function G_bd={&g1,11,0};

    size_t calls_bd=50000000;

    T_bd=gsl_rng_default;

    // bd, bulge to bulge
    double res2_bd,err2_bd;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2_bd[11]={3.0,3.0,0.03,-550.0,-450.0,-350.0,-550.0,-450.0,-350,0.0,-3.1415926/2.0};
    double xu2_bd[11]={13.0,13.0,0.08,550.0,450.0,350.0,550.0,450.0,350.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T2_bd;
    gsl_rng *r2_bd;

    gsl_monte_function G2_bd={&g2,11,0};

    size_t calls2_bd=50000000;

    T2_bd=gsl_rng_default;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // ms, disk to bulge
    double res_ms,err_ms;

    double xl_ms[11]={3.0,0.0,0.08,-100.0,-150.0,-100.0,-550.0,-450.0,-350,0.0,-3.1415926/2.0};
    double xu_ms[11]={13.0,13.0,1.0,100.0,150.0,100.0,550.0,450.0,350.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_ms;
    gsl_rng *r_ms;

    gsl_monte_function G_ms={&g1,11,0};

    size_t calls_ms=50000000;

    T_ms=gsl_rng_default;

    // ms, bulge to bulge
    double res2_ms,err2_ms;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2_ms[11]={3.0,3.0,0.08,-550.0,-450.0,-350.0,-550.0,-450.0,-350,0.0,-3.1415926/2.0};
    double xu2_ms[11]={13.0,13.0,1.0,550.0,450.0,350.0,550.0,450.0,350.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T2_ms;
    gsl_rng *r2_ms;

    gsl_monte_function G2_ms={&g2,11,0};

    size_t calls2_ms=50000000;

    T2_ms=gsl_rng_default;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // wd, disk to bulge
    double res_wd,err_wd;

    double xl_wd[11]={3.0,0.0,1.0,-100.0,-150.0,-100.0,-550.0,-450.0,-350,0.0,-3.1415926/2.0};
    double xu_wd[11]={13.0,13.0,8.0,100.0,150.0,100.0,550.0,450.0,350.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_wd;
    gsl_rng *r_wd;

    gsl_monte_function G_wd={&g1,11,0};

    size_t calls_wd=50000000;

    T_wd=gsl_rng_default;

    // wd, bulge to bulge
    double res2_wd,err2_wd;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2_wd[11]={3.0,3.0,1.0,-550.0,-450.0,-350.0,-550.0,-450.0,-350,0.0,-3.1415926/2.0};
    double xu2_wd[11]={13.0,13.0,8.0,550.0,450.0,350.0,550.0,450.0,350.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T2_wd;
    gsl_rng *r2_wd;

    gsl_monte_function G2_wd={&g2,11,0};

    size_t calls2_wd=50000000;

    T2_wd=gsl_rng_default;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // bh, disk to bulge
    double res_bh,err_bh;

    double xl_bh[11]={3.0,0.0,40.0,-100.0,-150.0,-100.0,-550.0,-450.0,-350,0.0,-3.1415926/2.0};
    double xu_bh[11]={13.0,13.0,120.0,100.0,150.0,100.0,550.0,450.0,350.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_bh;
    gsl_rng *r_bh;

    gsl_monte_function G_bh={&g1,11,0};

    size_t calls_bh=50000000;

    T_bh=gsl_rng_default;

    // bh, bulge to bulge
    double res2_bh,err2_bh;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2_bh[11]={3.0,3.0,40.0,-550.0,-450.0,-350.0,-550.0,-450.0,-350,0.0,-3.1415926/2.0};
    double xu2_bh[11]={13.0,13.0,120.0,550.0,450.0,350.0,550.0,450.0,350.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T2_bh;
    gsl_rng *r2_bh;

    gsl_monte_function G2_bh={&g2,11,0};

    size_t calls2_bh=50000000;

    T2_bh=gsl_rng_default;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // bd, disk to disk
    double res_bdd,err_bdd;

    double xl_bdd[11]={0.0,0.0,0.03,-100.0,-150.0,-100.0,-100.0,-150.0,-100,0.0,-3.1415926/2.0};
    double xu_bdd[11]={13.0,13.0,0.08,100.0,150.0,100.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_bdd;
    gsl_rng *r_bdd;

    gsl_monte_function G_bdd={&g1d,11,0};

    size_t calls_bdd=50000000;

    T_bdd=gsl_rng_default;

    // bd, bulge to disk
    double res2_bdd,err2_bdd;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2_bdd[11]={0.0,3.0,0.03,-550.0,-450.0,-350.0,-100.0,-150.0,-100,0.0,-3.1415926/2.0};
    double xu2_bdd[11]={13.0,13.0,0.08,550.0,450.0,350.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T2_bdd;
    gsl_rng *r2_bdd;

    gsl_monte_function G2_bdd={&g2d,11,0};

    size_t calls2_bdd=50000000;

    T2_bdd=gsl_rng_default;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // ms, disk to disk
    double res_msd,err_msd;

    double xl_msd[11]={0.0,0.0,0.08,-100.0,-150.0,-100.0,-100.0,-150.0,-100,0.0,-3.1415926/2.0};
    double xu_msd[11]={13.0,13.0,1.0,100.0,150.0,100.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_msd;
    gsl_rng *r_msd;

    gsl_monte_function G_msd={&g1d,11,0};

    size_t calls_msd=50000000;

    T_msd=gsl_rng_default;

    // ms, bulge to disk
    double res2_msd,err2_msd;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2_msd[11]={0.0,3.0,0.08,-550.0,-450.0,-350.0,-100.0,-150.0,-100,0.0,-3.1415926/2.0};
    double xu2_msd[11]={13.0,13.0,1.0,550.0,450.0,350.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T2_msd;
    gsl_rng *r2_msd;

    gsl_monte_function G2_msd={&g2d,11,0};

    size_t calls2_msd=50000000;

    T2_msd=gsl_rng_default;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // wd, disk to disk
    double res_wdd,err_wdd;

    double xl_wdd[11]={0.0,0.0,1.0,-100.0,-150.0,-100.0,-100.0,-150.0,-100,0.0,-3.1415926/2.0};
    double xu_wdd[11]={13.0,13.0,8.0,100.0,150.0,100.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_wdd;
    gsl_rng *r_wdd;

    gsl_monte_function G_wdd={&g1d,11,0};

    size_t calls_wdd=50000000;

    T_wdd=gsl_rng_default;

    // wd, bulge to disk
    double res2_wdd,err2_wdd;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2_wdd[11]={0.0,3.0,1.0,-550.0,-450.0,-350.0,-100.0,-150.0,-100,0.0,-3.1415926/2.0};
    double xu2_wdd[11]={13.0,13.0,8.0,550.0,450.0,350.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T2_wdd;
    gsl_rng *r2_wdd;

    gsl_monte_function G2_wdd={&g2d,11,0};

    size_t calls2_wdd=50000000;

    T2_wdd=gsl_rng_default;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // bh, disk to disk
    double res_bhd,err_bhd;

    double xl_bhd[11]={0.0,0.0,40.0,-100.0,-150.0,-100.0,-100.0,-150.0,-100,0.0,-3.1415926/2.0};
    double xu_bhd[11]={13.0,13.0,120.0,100.0,150.0,100.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_bhd;
    gsl_rng *r_bhd;

    gsl_monte_function G_bhd={&g1d,11,0};

    size_t calls_bhd=50000000;

    T_bhd=gsl_rng_default;

    // bh, bulge to disk
    double res2_bhd,err2_bhd;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2_bhd[11]={0.0,3.0,40.0,-550.0,-450.0,-350.0,-100.0,-150.0,-100,0.0,-3.1415926/2.0};
    double xu2_bhd[11]={13.0,13.0,120.0,550.0,450.0,350.0,100.0,150.0,100.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T2_bhd;
    gsl_rng *r2_bhd;

    gsl_monte_function G2_bhd={&g2d,11,0};

    size_t calls2_bhd=50000000;

    T2_bhd=gsl_rng_default;

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // source
    double res_s,err_s;

    double xl_s[3]={3.0,0.0,-3.1415926/2.0};
    double xu_s[3]={13.0,2.0*3.1415926,3.1415926/2.0};

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

    // source disk  
    double res_sd,err_sd;

    double xl_sd[3]={0.0,0.0,-3.1415926/2.0};
    double xu_sd[3]={13.0,2.0*3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_sd;
    gsl_rng *r_sd;

    gsl_monte_function G_sd={&source_disk,3,0};

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

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G, xl, xu, 11, 10000, r, s,&res, &err);
    //display_results ("vegas warm-up", res, err);
    //printf ("converging...\n");

    do
    {
        gsl_monte_vegas_integrate (&G, xl, xu, 11, calls/5, r, s,&res, &err);
        //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    //display_results ("vegas final", res, err);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r2 = gsl_rng_alloc (T2);

    gsl_monte_vegas_state *s2 = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G2, xl2, xu2, 11, 10000, r2, s2, &res2, &err2);

    do
    {
        gsl_monte_vegas_integrate (&G2, xl2, xu2, 11, calls2/5, r2, s2, &res2, &err2);
    }
    while (fabs (gsl_monte_vegas_chisq (s2) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r_bd = gsl_rng_alloc (T_bd);

    gsl_monte_vegas_state *s_bd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G_bd, xl_bd, xu_bd, 11, 10000, r_bd, s_bd, &res_bd, &err_bd);

    do
    {
        gsl_monte_vegas_integrate (&G_bd, xl_bd, xu_bd, 11, calls_bd/5, r_bd, s_bd, &res_bd, &err_bd);
    }
    while (fabs (gsl_monte_vegas_chisq (s_bd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r2_bd = gsl_rng_alloc (T2_bd);

    gsl_monte_vegas_state *s2_bd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G2_bd, xl2_bd, xu2_bd, 11, 10000, r2_bd, s2_bd, &res2_bd, &err2_bd);

    do
    {
        gsl_monte_vegas_integrate (&G2_bd, xl2_bd, xu2_bd, 11, calls2_bd/5, r2_bd, s2_bd, &res2_bd, &err2_bd);
    }
    while (fabs (gsl_monte_vegas_chisq (s2_bd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r_ms = gsl_rng_alloc (T_ms);

    gsl_monte_vegas_state *s_ms = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G_ms, xl_ms, xu_ms, 11, 10000, r_ms, s_ms, &res_ms, &err_ms);

    do
    {
        gsl_monte_vegas_integrate (&G_ms, xl_ms, xu_ms, 11, calls_ms/5, r_ms, s_ms, &res_ms, &err_ms);
    }
    while (fabs (gsl_monte_vegas_chisq (s_ms) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r2_ms = gsl_rng_alloc (T2_ms);

    gsl_monte_vegas_state *s2_ms = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G2_ms, xl2_ms, xu2_ms, 11, 10000, r2_ms, s2_ms, &res2_ms, &err2_ms);

    do
    {
        gsl_monte_vegas_integrate (&G2_ms, xl2_ms, xu2_ms, 11, calls2_ms/5, r2_ms, s2_ms, &res2_ms, &err2_ms);
    }
    while (fabs (gsl_monte_vegas_chisq (s2_ms) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r_wd = gsl_rng_alloc (T_wd);

    gsl_monte_vegas_state *s_wd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G_wd, xl_wd, xu_wd, 11, 10000, r_wd, s_wd, &res_wd, &err_wd);

    do
    {
        gsl_monte_vegas_integrate (&G_wd, xl_wd, xu_wd, 11, calls_wd/5, r_wd, s_wd, &res_wd, &err_wd);
    }
    while (fabs (gsl_monte_vegas_chisq (s_wd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r2_wd = gsl_rng_alloc (T2_wd);

    gsl_monte_vegas_state *s2_wd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G2_wd, xl2_wd, xu2_wd, 11, 10000, r2_wd, s2_wd, &res2_wd, &err2_wd);

    do
    {
        gsl_monte_vegas_integrate (&G2_wd, xl2_wd, xu2_wd, 11, calls2_wd/5, r2_wd, s2_wd, &res2_wd, &err2_wd);
    }
    while (fabs (gsl_monte_vegas_chisq (s2_wd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r_bh = gsl_rng_alloc (T_bh);

    gsl_monte_vegas_state *s_bh = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G_bh, xl_bh, xu_bh, 11, 10000, r_bh, s_bh, &res_bh, &err_bh);

    do
    {
        gsl_monte_vegas_integrate (&G_bh, xl_bh, xu_bh, 11, calls_bh/5, r_bh, s_bh, &res_bh, &err_bh);
    }
    while (fabs (gsl_monte_vegas_chisq (s_bh) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r2_bh = gsl_rng_alloc (T2_bh);

    gsl_monte_vegas_state *s2_bh = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G2_bh, xl2_bh, xu2_bh, 11, 10000, r2_bh, s2_bh, &res2_bh, &err2_bh);

    do
    {
        gsl_monte_vegas_integrate (&G2_bh, xl2_bh, xu2_bh, 11, calls2_bh/5, r2_bh, s2_bh, &res2_bh, &err2_bh);
    }
    while (fabs (gsl_monte_vegas_chisq (s2_bh) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    rd = gsl_rng_alloc (Td);

    gsl_monte_vegas_state *sd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&Gd, xld, xud, 11, 10000, rd, sd,&resd, &errd);
    //display_results ("vegas warm-up", res, err);
    //printf ("converging...\n");

    do
    {
        gsl_monte_vegas_integrate (&Gd, xld, xud, 11, callsd/5, rd, sd,&resd, &errd);
        //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
    }
    while (fabs (gsl_monte_vegas_chisq (sd) - 1.0) > 0.5);
    //display_results ("vegas final", res, err);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r2d = gsl_rng_alloc (T2d);

    gsl_monte_vegas_state *s2d = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G2d, xl2d, xu2d, 11, 10000, r2d, s2d, &res2d, &err2d);

    do
    {
        gsl_monte_vegas_integrate (&G2d, xl2d, xu2d, 11, calls2d/5, r2d, s2d, &res2d, &err2d);
    }
    while (fabs (gsl_monte_vegas_chisq (s2d) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r_bdd = gsl_rng_alloc (T_bdd);

    gsl_monte_vegas_state *s_bdd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G_bdd, xl_bdd, xu_bdd, 11, 10000, r_bdd, s_bdd, &res_bdd, &err_bdd);

    do
    {
        gsl_monte_vegas_integrate (&G_bdd, xl_bdd, xu_bdd, 11, calls_bdd/5, r_bdd, s_bdd, &res_bdd, &err_bdd);
    }
    while (fabs (gsl_monte_vegas_chisq (s_bdd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r2_bdd = gsl_rng_alloc (T2_bdd);

    gsl_monte_vegas_state *s2_bdd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G2_bdd, xl2_bdd, xu2_bdd, 11, 10000, r2_bdd, s2_bdd, &res2_bdd, &err2_bdd);

    do
    {
        gsl_monte_vegas_integrate (&G2_bdd, xl2_bdd, xu2_bdd, 11, calls2_bdd/5, r2_bdd, s2_bdd, &res2_bdd, &err2_bdd);
    }
    while (fabs (gsl_monte_vegas_chisq (s2_bdd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r_msd = gsl_rng_alloc (T_msd);

    gsl_monte_vegas_state *s_msd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G_msd, xl_msd, xu_msd, 11, 10000, r_msd, s_msd, &res_msd, &err_msd);

    do
    {
        gsl_monte_vegas_integrate (&G_msd, xl_msd, xu_msd, 11, calls_msd/5, r_msd, s_msd, &res_msd, &err_msd);
    }
    while (fabs (gsl_monte_vegas_chisq (s_msd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r2_msd = gsl_rng_alloc (T2_msd);

    gsl_monte_vegas_state *s2_msd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G2_msd, xl2_msd, xu2_msd, 11, 10000, r2_msd, s2_msd, &res2_msd, &err2_msd);

    do
    {
        gsl_monte_vegas_integrate (&G2_msd, xl2_msd, xu2_msd, 11, calls2_msd/5, r2_msd, s2_msd, &res2_msd, &err2_msd);
    }
    while (fabs (gsl_monte_vegas_chisq (s2_msd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r_wdd = gsl_rng_alloc (T_wdd);

    gsl_monte_vegas_state *s_wdd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G_wdd, xl_wdd, xu_wdd, 11, 10000, r_wdd, s_wdd, &res_wdd, &err_wdd);

    do
    {
        gsl_monte_vegas_integrate (&G_wdd, xl_wdd, xu_wdd, 11, calls_wdd/5, r_wdd, s_wdd, &res_wdd, &err_wdd);
    }
    while (fabs (gsl_monte_vegas_chisq (s_wdd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r2_wdd = gsl_rng_alloc (T2_wdd);

    gsl_monte_vegas_state *s2_wdd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G2_wdd, xl2_wdd, xu2_wdd, 11, 10000, r2_wdd, s2_wdd, &res2_wdd, &err2_wdd);

    do
    {
        gsl_monte_vegas_integrate (&G2_wdd, xl2_wdd, xu2_wdd, 11, calls2_wdd/5, r2_wdd, s2_wdd, &res2_wdd, &err2_wdd);
    }
    while (fabs (gsl_monte_vegas_chisq (s2_wdd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r_bhd = gsl_rng_alloc (T_bhd);

    gsl_monte_vegas_state *s_bhd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G_bhd, xl_bhd, xu_bhd, 11, 10000, r_bhd, s_bhd, &res_bhd, &err_bhd);

    do
    {
        gsl_monte_vegas_integrate (&G_bhd, xl_bhd, xu_bhd, 11, calls_bhd/5, r_bhd, s_bhd, &res_bhd, &err_bhd);
    }
    while (fabs (gsl_monte_vegas_chisq (s_bhd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r2_bhd = gsl_rng_alloc (T2_bhd);

    gsl_monte_vegas_state *s2_bhd = gsl_monte_vegas_alloc (11);
    gsl_monte_vegas_integrate (&G2_bhd, xl2_bhd, xu2_bhd, 11, 10000, r2_bhd, s2_bhd, &res2_bhd, &err2_bhd);

    do
    {
        gsl_monte_vegas_integrate (&G2_bhd, xl2_bhd, xu2_bhd, 11, calls2_bhd/5, r2_bhd, s2_bhd, &res2_bhd, &err2_bhd);
    }
    while (fabs (gsl_monte_vegas_chisq (s2_bhd) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
    
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

    r_xd=gsl_rng_alloc(T_xd);
        
    gsl_monte_vegas_state *s_xd = gsl_monte_vegas_alloc (10);
    gsl_monte_vegas_integrate (&G_xd, xl_xd, xu_xd, 10, 10000, r_xd, s_xd, &res_xd, &err_xd);

    do
    {
        gsl_monte_vegas_integrate (&G_xd, xl_xd, xu_xd, 10, calls_xd/5, r_xd, s_xd, &res_xd, &err_xd);
    }
    while (fabs (gsl_monte_vegas_chisq (s_xd) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    bd = (fabs(res_bd)+fabs(res_bdd)+fabs(res2_bd)+fabs(res2_bdd))/(fabs(res_s)+fabs(res_sd));
    ms = (fabs(res_ms)+fabs(res_msd)+fabs(res2_ms)+fabs(res2_msd))/(fabs(res_s)+fabs(res_sd));
    wd = (fabs(res_wd)+fabs(res_wdd)+fabs(res2_wd)+fabs(res2_wdd))/(fabs(res_s)+fabs(res_sd));
    bh = (fabs(res_bh)+fabs(res_bhd)+fabs(res2_bh)+fabs(res2_bhd))/(fabs(res_s)+fabs(res_sd));
    ns = (fabs(res_x)+fabs(res_xd))/(fabs(res_s)+fabs(res_sd));
    total = (fabs(res)+fabs(resd)+fabs(res2)+fabs(res2d))/(fabs(res_s)+fabs(res_sd))+ns;
    ns_p = ns/total;
    bd_p = bd/total;
    ms_p = ms/total;
    wd_p = wd/total;
    bh_p = bh/total;
    printf ("bd: %e\n", bd_p);
    printf ("ms: %e\n", ms_p);
    printf ("wd: %e\n", wd_p);
    printf ("bh: %e\n", bh_p);
    printf ("ns: %e\n", ns_p);

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    gsl_monte_vegas_free (s2);

    gsl_rng_free (r2);

    gsl_monte_vegas_free (s_bd);

    gsl_rng_free (r_bd);

    gsl_monte_vegas_free (s2_bd);

    gsl_rng_free (r2_bd);

    gsl_monte_vegas_free (s_ms);

    gsl_rng_free (r_ms);

    gsl_monte_vegas_free (s2_ms);

    gsl_rng_free (r2_ms);

    gsl_monte_vegas_free (s_wd);

    gsl_rng_free (r_wd);

    gsl_monte_vegas_free (s2_wd);

    gsl_rng_free (r2_wd);

    gsl_monte_vegas_free (s_bh);

    gsl_rng_free (r_bh);

    gsl_monte_vegas_free (s2_bh);

    gsl_rng_free (r2_bh);

    gsl_monte_vegas_free (sd);

    gsl_rng_free (rd);

    gsl_monte_vegas_free (s2d);

    gsl_rng_free (r2d);

    gsl_monte_vegas_free (s_bdd);

    gsl_rng_free (r_bdd);

    gsl_monte_vegas_free (s2_bdd);

    gsl_rng_free (r2_bdd);

    gsl_monte_vegas_free (s_msd);

    gsl_rng_free (r_msd);

    gsl_monte_vegas_free (s2_msd);

    gsl_rng_free (r2_msd);

    gsl_monte_vegas_free (s_wdd);

    gsl_rng_free (r_wdd);

    gsl_monte_vegas_free (s2_wdd);

    gsl_rng_free (r2_wdd);

    gsl_monte_vegas_free (s_bhd);

    gsl_rng_free (r_bhd);

    gsl_monte_vegas_free (s2_bhd);

    gsl_rng_free (r2_bhd);

    gsl_monte_vegas_free (s_x);

    gsl_rng_free (r_x);

    gsl_monte_vegas_free (s_xd);

    gsl_rng_free (r_xd);

    return 0;
}
