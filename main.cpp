#include <iostream>
#include <fstream>
#include <iomanip>
#include "RK4Solver.hpp"

int main() {

    // define initial data and stepsize 
    double r0 = 0.0; // start radius 
    double r1 = 1000.0; // end radius 
    double lapse_0 = 1.0;  // Initial lapse
    double psi_0 = 1.0;  // Initial conformal factor
    double dr = 0.004; // stepsize

    // define the integrator
    RK4Solver solver;

    // result, vector of tuples
    auto result = solver.integrate(r0, r1, lapse_0, psi_0, dr);

    // find the asymptotic value of the lapse
    // then re-scale the initial guess lapse_0
    lapse_0 = lapse_0/solver.renorm;

    // re-integrate with better aymptotics
    result = solver.integrate(r0, r1, lapse_0, psi_0, dr);

    // Open output file
    std::ofstream outFile("output.csv");
    if (!outFile) {
        std::cerr << "Error: Unable to open file for writing.\n";
        return 1;
    }

    // Write header
   outFile << "r," << "alpha_norm(r)," << "psi(r)," << "u(r)\n";

    // Write data
    outFile << std::fixed << std::setprecision(5);
    for (const auto& [r, lapse, psi, u] : result) {
        outFile << r << "," << lapse << "," << psi << "," << u << "\n";
    }

    outFile.close();
    std::cout << "Data written to output.csv\n";


    std::cout << "Total Mass : " << solver.mass_inf << std::endl;

    return 0;
}
