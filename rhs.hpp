#ifndef RHS_HPP
#define RHS_HPP

#include <cmath>

class RHS {
public:
    // constructor 
    RHS();
    
    // gives the scalar field output
    double lapse(double r, double lapse, double psi, double u);
    double psi(double r, double lapse, double psi, double u);

private:
    double eps = pow(10.,-10);
};

#endif // RHS_HPP
