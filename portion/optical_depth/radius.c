// calculate the percentage contribution of optical depth from different types of lens 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main (void)
{
    double F,ita,r;

    // lense distribution: disk

    double p0=0.01*1.3888*0.0493*pow(10,9); // M*kpc^(-3)
    double beta=0.565;
    double h1=0.270; // kpc
    double h2=0.440; // kpc
    double H=2.75; // kpc
    double z=0.0;
    double r0=8.0;

    for (r=0.0;r<=20.0;r=r+0.01)
    {
        if (((r/9.025)+0.114)<=0.670)
        {
            ita=0.670;
        }
        else
        {
            ita=(r/9.025)+0.114;
        }

        F=(p0/ita)*exp(-(r-r0)/H)*((1.0-beta)*pow(cosh(z/(ita*h1)),-2.0)+beta*exp(-fabs(z)/(ita*h2)));  // :disk
        printf("%e %e\n", r, F);
    }

    return 0;
}

