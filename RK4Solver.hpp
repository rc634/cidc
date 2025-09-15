#ifndef RK4SOLVER_HPP
#define RK4SOLVER_HPP

#include <functional>
#include <vector>
#include <utility>
#include "scalar_profile.hpp"
#include "rhs.hpp"

class RK4Solver {
public:
    //constructor
    RK4Solver();

    // Integrate from r0 to r1 with step size h
    std::vector<std::tuple<double, double, double, double, double>> 
    integrate(  double r0, double r1,
                double lapse0, double psi0, double zeta0,
                double h);

    // Adaptive integrate from r0 to r2 with step size h, 
    //adaptive stepsize triggers at 'knee'
    std::vector<std::tuple<double, double, double, double, double>> 
    integrate_adaptive(  double r0, double knee, double r2,
                         double lapse0, double psi0, double zeta0,
                         double h);

private:
    Scalar_Profile scalar;
    RHS rhs;

public:
    // large radius values
    double lapse_inf = 0.;
    double psi_inf = 0.;
    double zeta_inf = 0.;
    // mass of solution
    double mass_inf = 0.; 
    // used to renormalise alpha
    double renorm = 0.;
};

#endif // RK4SOLVER_HPP
