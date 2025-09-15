#include "RK4Solver.hpp"
#include <iostream>

RK4Solver::RK4Solver()
{
    // empty const
}

std::vector<std::tuple<double, double, double, double, double>> 
RK4Solver::integrate(
    double r0, double r1,
    double lapse0, double psi0, double zeta0,
    double h)
{
    std::vector<std::tuple<double, double, double, double, double>> result;
    double r = r0;
    double lapse = lapse0;
    double psi = psi0;
    double zeta = zeta0;

    result.emplace_back(r, lapse, psi, zeta, scalar.u_of_r(r));

    while (r < r1) {
        // Adjust last step to not overshoot
        if (r + h > r1) 
        {
            h = r1 - r;
            lapse_inf = lapse;
            psi_inf = psi;
            mass_inf = 2. * r * (psi_inf - 1.);
        }

        double up1 = scalar.du_dr(r);
        double k1 = h * rhs.lapse(r, lapse, psi, zeta, up1);
        double q1 = h * rhs.psi(r, lapse, psi, zeta, up1);
        double z1 = h * rhs.zeta(r, lapse, psi, zeta, up1);

        double up2 = scalar.du_dr(r + h/2);
        double k2 = h * rhs.lapse(r + h/2, lapse + k1/2, psi + q1/2, zeta + z1/2, up2);
        double q2 = h * rhs.psi(r + h/2, lapse + k1/2, psi + q1/2, zeta + z1/2, up2);
        double z2 = h * rhs.zeta(r + h/2, lapse + k1/2, psi + q1/2, zeta + z1/2, up2);

        double up3 = scalar.du_dr(r + h/2);
        double k3 = h * rhs.lapse(r + h/2, lapse + k2/2, psi + q2/2, zeta + z2/2, up3);
        double q3 = h * rhs.psi(r + h/2, lapse + k2/2, psi + q2/2, zeta + z2/2, up3);
        double z3 = h * rhs.zeta(r + h/2, lapse + k2/2, psi + q2/2, zeta + z2/2, up3);

        double up4 = scalar.du_dr(r + h);
        double k4 = h * rhs.lapse(r + h, lapse + k3, psi + q3, zeta + z3, up4);
        double q4 = h * rhs.psi(r + h, lapse + k3, psi + q3, zeta + z3, up4);
        double z4 = h * rhs.zeta(r + h, lapse + k3, psi + q3, zeta + z3, up4);

        lapse += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
        psi += (q1 + 2*q2 + 2*q3 + q4) / 6.0;
        zeta += (z1 + 2*z2 + 2*z3 + z4) / 6.0;
        r += h;

        result.emplace_back(r, lapse, psi, zeta, scalar.u_of_r(r));
    }

    return result;
}

std::vector<std::tuple<double, double, double, double, double>> 
RK4Solver::integrate_adaptive(
    double r0, double knee, double r1,
    double lapse0, double psi0, double zeta0,
    double h)
{
    std::vector<std::tuple<double, double, double, double, double>> result;
    double r = r0;
    double lapse = lapse0;
    double psi = psi0;
    double zeta = zeta0;
    double h0 = h; // initial h or dr saved for adaptive part

    result.emplace_back(r, lapse, psi, zeta, scalar.u_of_r(r));

    while (r < r1) {
        // Adjust last step to not overshoot
        if (r + h > r1) 
        {
            h = r1 - r;
            lapse_inf = lapse;
            psi_inf = psi;
            mass_inf = 2. * r * (psi_inf - 1.);
        }

        if (r>knee)
        {
            // increases stepsize after 500 - hard coded here
            h = r/knee * h0;
        }

        double up1 = scalar.du_dr(r);
        double k1 = h * rhs.lapse(r, lapse, psi, zeta, up1);
        double q1 = h * rhs.psi(r, lapse, psi, zeta, up1);
        double z1 = h * rhs.zeta(r, lapse, psi, zeta, up1);

        double up2 = scalar.du_dr(r + h/2);
        double k2 = h * rhs.lapse(r + h/2, lapse + k1/2, psi + q1/2, zeta + z1/2, up2);
        double q2 = h * rhs.psi(r + h/2, lapse + k1/2, psi + q1/2, zeta + z1/2, up2);
        double z2 = h * rhs.zeta(r + h/2, lapse + k1/2, psi + q1/2, zeta + z1/2, up2);

        double up3 = scalar.du_dr(r + h/2);
        double k3 = h * rhs.lapse(r + h/2, lapse + k2/2, psi + q2/2, zeta + z2/2, up3);
        double q3 = h * rhs.psi(r + h/2, lapse + k2/2, psi + q2/2, zeta + z2/2, up3);
        double z3 = h * rhs.zeta(r + h/2, lapse + k2/2, psi + q2/2, zeta + z2/2, up3);

        double up4 = scalar.du_dr(r + h);
        double k4 = h * rhs.lapse(r + h, lapse + k3, psi + q3, zeta + z3, up4);
        double q4 = h * rhs.psi(r + h, lapse + k3, psi + q3, zeta + z3, up4);
        double z4 = h * rhs.zeta(r + h, lapse + k3, psi + q3, zeta + z3, up4);

        lapse += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
        psi += (q1 + 2*q2 + 2*q3 + q4) / 6.0;
        zeta += (z1 + 2*z2 + 2*z3 + z4) / 6.0;
        r += h;

        result.emplace_back(r, lapse, psi, zeta, scalar.u_of_r(r));
    }

    return result;
}
