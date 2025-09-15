#ifndef RHS_HPP
#define RHS_HPP

#include <cmath>

class RHS {
public:
    // constructor 
    RHS();
    
    // gives the scalar field output
    double lapse(double r, double lapse, double psi, double zeta, double up);
    double psi(double r, double lapse, double psi, double zeta, double up);
    double zeta(double r, double lapse, double psi, double zeta, double up);


private:
    double eps = pow(10.,-10);
};

#endif // RHS_HPP
