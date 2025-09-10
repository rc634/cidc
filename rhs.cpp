
#include "rhs.hpp"
#include <cmath>

// empty constructor 
RHS::RHS()
{
    // empty constructor
}


// gerneal note
// up = u prime = u' = du/dr


// return scalar field 
double RHS::lapse(double r, double lapse, double psi, double up)
{
    double safe_r = r;
    if (r < eps)
    {
        safe_r = eps;
    }
    double lapse_rhs = lapse*(-1. + psi*psi + 4.*M_PI*r*r*up*up)/(2.*safe_r);
    return lapse_rhs;
}

// return scalar field 
double RHS::psi(double r, double lapse, double psi, double up)
{
    double safe_r = r;
    if (r < eps)
    {
        safe_r = eps;
    }
    double psi_rhs = psi*(1. - psi*psi + 4.*M_PI*r*r*up*up)/(2.*safe_r);
    return psi_rhs;
}
