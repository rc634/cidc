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
    double gaussian = u0 * exp(-x*x/sigma);
    double u = a_r * a_r * a_r * gaussian;
    return u;
}

// return scalar field deriv
double Scalar_Profile::du_dr(double a_r)
{
    double x = a_r - rc;
    double gaussian = u0 * exp(-x*x/sigma);
    double dudr = (3. - 2. * x * a_r/sigma)*a_r*a_r*gaussian;
    return dudr;
}
