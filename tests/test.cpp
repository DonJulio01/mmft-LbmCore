#include <gtest/gtest.h>
#include "poiseuille_helper.h"
#include <iostream>

// Tests if the numerical solution is close enough to the analytical solution
TEST(PoiseuilleTest, ProfileAccuracy) {
    std::size_t nx    = 200;
    std::size_t ny    = 100;
    double      tau   = 0.7;
    double      u_max = 0.05;
    std::size_t steps = 50000;
    double expected_error = 1e-3;

    double error = compute_poiseuille_error(nx, ny, tau, u_max, steps);

    std::cout << "Error: " << error << std::endl;
    std::cout << "Expected error: " << expected_error << std::endl;

    // Allowed tolerance
    EXPECT_LT(error, expected_error);
}
