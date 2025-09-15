#ifndef SCALAR_PROFILE_HPP
#define SCALAR_PROFILE_HPP

class Scalar_Profile {
public:
    // constructor 
    Scalar_Profile();
    
    // gives the scalar field output
    double u_of_r(double a_r);

    // gives the scalar field derivative
    double du_dr(double a_r);

    // gives the scalar field 2nd derivative
    double d2u_dr2(double a_r);

private:
    double u0 = 0.00005;
    double sigma = 1.;
    double rc = 10.;
};

#endif // SCALAR_PROFILE_HPP
