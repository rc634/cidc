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
    std::vector<std::tuple<double, double, double, double>> 
    integrate(  double r0, double r1,
                double lapse0, double psi0,
                double h);

private:
    Scalar_Profile scalar;
    RHS rhs;

public:
    // large radius value of alpha
    double alpha_inf = 1.;
};

#endif // RK4SOLVER_HPP
