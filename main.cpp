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

    // result tuple 
    auto result = solver.integrate(r0, r1, lapse_0, psi_0, dr);

    // print to screen 
    // std::cout << std::fixed << std::setprecision(5);
    // std::cout << "Time\tY1\t\tY2\n";
    // for (const auto& [t, y1, y2] : result) {
    //     std::cout << t << "\t" << y1 << "\t" << y2 << "\n";
    // }

    // find the final value of alpha
    lapse_0 = lapse_0/solver.alpha_inf;
    
    // re-integrate with better aymptotics
    result = solver.integrate(r0, r1, lapse_0, psi_0, dr);

    // Open output file
    std::ofstream outFile("output.csv");
    if (!outFile) {
        std::cerr << "Error: Unable to open file for writing.\n";
        return 1;
    }

    // Write header
    outFile << "Radius," << "Lapse," << "Psi," << "U\n";

    // Write data
    outFile << std::fixed << std::setprecision(5);
    for (const auto& [r, lapse, psi, u] : result) {
        outFile << r << "," << lapse << "," << psi << "," << u << "\n";
    }

    outFile.close();
    std::cout << "Data written to output.csv\n";

    return 0;
}
