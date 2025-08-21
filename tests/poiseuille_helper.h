#pragma once

#include "../src/lbm.h"

// Returns error between numerical and analytical solution
inline double compute_poiseuille_error(std::size_t nx,
                                       std::size_t ny,
                                       double tau,
                                       double u_max,
                                       std::size_t steps) {

    lbm::Solver solver(nx, ny, tau, u_max, steps);
    solver.run();
    return solver.computeError();
}
