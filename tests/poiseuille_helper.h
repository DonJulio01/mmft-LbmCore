#pragma once

#include "../src/lbm.h"

// Returns error between numerical and analytical solution
inline double compute_poiseuille_error(std::size_t nx,
                                       std::size_t ny,
                                       double tau,
                                       double u_max,
                                       std::size_t steps) {

    double nu = (tau - 0.5) / 3.0;
    double H  = static_cast<double>(ny - 2);
    double force_x = 0.0;
    if (H > 0.0) {
        force_x = 8.0 * u_max * nu / (H * H);
    }

    lbm::Solver solver(nx, ny, tau, force_x, u_max, steps);
    solver.run();
    return solver.computeError();
}
