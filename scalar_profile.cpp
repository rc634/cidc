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

// return scalar field deriv
double Scalar_Profile::d2u_dr2(double a_r)
{
    double x = a_r - rc;
    double gaussian = u0 * exp(-x*x/sigma);
    double sigsqr = sigma*sigma;
    double d2udr2 = 6.*a_r + 4.*pow(a_r,5)/sigsqr  
                  - 8.*pow(a_r,4)*rc/sigsqr + 4.*pow(a_r,3)*rc*rc/sigsqr 
                  - 14.*pow(a_r,3)/sigma + 12*a_r*a_r*rc/sigma;
    d2udr2 *= a_r*gaussian;
    return d2udr2;
}

