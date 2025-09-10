#include "scalar_profile.hpp"
#include <cmath>

// empty constructor 
Scalar_Profile::Scalar_Profile()
{
    // empty constructor 
}

// return scalar field 
double Scalar_Profile::u_of_r(double a_r)
{
    double x = a_r - rc;
    double u = u0 * exp(-x*x/sigma);
    return u;
}

// return scalar field deriv
double Scalar_Profile::du_dr(double a_r)
{
    double x = a_r - rc;
    double dudr = (- 2. * x /sigma) * u_of_r(a_r);
    return dudr;
}
