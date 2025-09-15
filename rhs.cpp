
#include "rhs.hpp"
#include <cmath>

// empty constructor 
RHS::RHS()
{
    // empty constructor
}


// gerneal note
// up = u prime = u' = du/dr


// calculate RHS for lapse
double RHS::lapse(double r, double lapse, double psi, double zeta, double up)
{
    double safe_r = r;
    if (r < eps)
    {
        safe_r = eps;
    }
    double numerator = M_PI*r*psi*psi*up*up - psi*zeta - r*zeta*zeta;
    double denominator = psi*(psi + 2. * r * zeta);
    return 2. * lapse * numerator / denominator;
}

// calculate RHS for psi
double RHS::psi(double r, double lapse, double psi, double zeta, double up)
{
    double safe_r = r;
    if (r < eps)
    {
        safe_r = eps;
    }
    double psi_rhs = zeta;
    return psi_rhs;
}

// calculate RHS for zeta or psi'
double RHS::zeta(double r, double lapse, double psi, double zeta, double up)
{
    double safe_r = r;
    if (r < eps)
    {
        safe_r = eps;
    }
    double zeta_rhs = - M_PI * psi * up * up - 2. * zeta / safe_r;
    return zeta_rhs;
}
