#include <iostream>
#include <cstdlib>

#include "lbm.h"

int main(int argc, char const* argv[]) {
    std::size_t nx = 100;
    std::size_t ny = 50;
    double tau = 0.7;
    double u_max = 0.05;
    std::size_t steps = 10000;
    if (argc >= 2) {
        nx = static_cast<std::size_t>(std::atoi(argv[1]));
    }
    if (argc >= 3) {
        ny = static_cast<std::size_t>(std::atoi(argv[2]));
    }
    if (argc >= 4) {
        tau = std::atof(argv[3]);
    }
    if (argc >= 5) {
        u_max = std::atof(argv[4]);
    }
    if (argc >= 6) {
        steps = static_cast<std::size_t>(std::atoll(argv[5]));
    }

    // Compute the body force needed to achieve the desired maximum velocity for plane Poiseuille flow.
    double nu = (tau - 0.5) / 3.0;
    double H = static_cast<double>(ny - 2);
    double force_x = 0.0;
    if (H > 0.0) {
        // From analytical steady-state Poiseuille formula
        force_x = 8.0 * u_max * nu / (H * H);
    }

    std::cout << "Running Poiseuille flow simulation on a " << nx << "x" << ny
              << " lattice for " << steps << " time steps..." << std::endl;
    std::cout << "tau = " << tau << ", u_max = " << u_max
              << ", computed force_x = " << force_x << std::endl;

    try {
        lbm::Solver solver(nx, ny, tau, force_x, u_max, steps);
        solver.run();
        // Write the resulting velocity profile to disk
        solver.writeVelocityProfile("velocity_profile.csv");
        std::cout << "Resulting velocity profile successfully written to disk." << std::endl;
    } catch (const std::exception &ex) {
        std::cerr << "Error during simulation: " << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
    return 0;
}
