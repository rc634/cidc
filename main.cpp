#include <iostream>
#include <fstream>
#include <iomanip>
#include "RK4Solver.hpp"

int main() {

    // define initial data and stepsize 
    double r0 = 0.0; // start radius 
    double r1 = 2000.0; // end radius -- used as knee in adaptive integrator too
    double r2 = 2000000.0; // asymptotics radius (numerical infinity)
    double lapse_0 = 1.0;  // Initial lapse
    double psi_0 = 1.0;  // Initial conformal factor
    double zeta_0 = 0.0; // Initial conformal factor derivative - MUST be zero
    double dr = 0.004; // stepsize

    // define the integrator
    RK4Solver solver;

    // result, vector of tuples
    auto result = solver.integrate_adaptive(r0, r1, r2, lapse_0, psi_0, zeta_0, dr);

    // find the asymptotic value of the lapse
    // then re-scale the initial guess lapse_0
    lapse_0 = lapse_0/solver.lapse_inf;
    psi_0 = psi_0/solver.psi_inf;


    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "-- * STARTING CHOPTUIK INITIAL DATA SOLVER * --" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "-- DO ASYMPTOTIC RK4" << std::endl;


    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "Adaptive rk4 to r = " << r2 << " found asymptotics :" << std::endl;
    std::cout << "Psi_inf = " << solver.psi_inf << " and lapse_inf = " << solver.lapse_inf << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    // re-integrate with better aymptotics
    result = solver.integrate(r0, r1, lapse_0, psi_0, zeta_0, dr);

    std::cout << "-- DO REGULAR RK4" << std::endl;

    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "Normal rk4 to r = " << r1 << " found asymptotics :" << std::endl;
    std::cout << "Psi_inf = " << solver.psi_inf << " and lapse_inf = " << solver.lapse_inf << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    // Open output file
    std::ofstream outFile("output.csv");
    if (!outFile) {
        std::cerr << "Error: Unable to open file for writing.\n";
        return 1;
    }

    // Write header
   outFile << "r," << "alpha_norm(r)," << "psi(r)," << "zeta(r)," << "u(r)\n";

    // Write data
    outFile << std::fixed << std::setprecision(5);
    for (const auto& [r, lapse, psi, zeta, u] : result) {
        outFile << r << "," << lapse << "," << psi << "," << zeta << "," << u << "\n";
    }

    outFile.close();
    std::cout << "Data written to output.csv\n";


    std::cout << "Total Mass : " << solver.mass_inf << std::endl;

    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "-- * END * --" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    return 0;
}
