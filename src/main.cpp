#include <iostream>
#include <cstdlib>

#include "lbm.h"

/**
 * @param argv
 *   argv[1] = w_phys [m]           — channel width
 *   argv[2] = nu_phys [m^2/s]      — kinematic viscosity
 *   argv[3] = rho_phys [kg/m^3]    — density
 *   argv[4] = u_max_phys [m/s]     — target peak velocity
 *   argv[5] = u_max_lu             — target peak velocity in lattice units (for numerical stability)
 *   argv[6] = nx [cells]           — x-direction cells
 *   argv[7] = ny [cells]           — y-direction cells (incl. walls)
 *   argv[8] = steps                — time steps (physical duration = steps * dt)
 */
int main(int argc, char const* argv[]) {
    if (argc < 9) {
        std::cerr << "usage: w[m] nu[m2/s] rho[kg/m3] u_max[m/s] u_max_lu nx ny steps\n";
        return EXIT_FAILURE;
    }
    const double w_phys   = std::atof(argv[1]);
    const double nu_phys  = std::atof(argv[2]);
    const double rho_phys = std::atof(argv[3]);
    const double u_max_p  = std::atof(argv[4]);
    const double u_max_lu = std::atof(argv[5]);
    const std::size_t nx  = static_cast<std::size_t>(std::atoi(argv[6]));
    const std::size_t ny  = static_cast<std::size_t>(std::atoi(argv[7]));
    const std::size_t steps = static_cast<std::size_t>(std::atoi(argv[8]));

    if (!(w_phys  > 0.0)) { std::cerr << "w_phys must be > 0\n"; return EXIT_FAILURE; }
    if (!(nu_phys > 0.0)) { std::cerr << "nu_phys must be > 0\n"; return EXIT_FAILURE; }
    if (!(rho_phys> 0.0)) { std::cerr << "rho_phys must be > 0\n"; return EXIT_FAILURE; }
    if (!(u_max_p > 0.0)) { std::cerr << "u_max_phys must be > 0\n"; return EXIT_FAILURE; }
    if (!(u_max_lu > 0.0 && u_max_lu < 0.2)) { std::cerr << "u_max_lu must be in (0,0.2)\n"; return EXIT_FAILURE; }
    if (nx < 8) { std::cerr << "nx too small (<8)\n"; return EXIT_FAILURE; }
    if (ny <= 2) { std::cerr << "ny must be > 2\n"; return EXIT_FAILURE; }
    if (!(steps > 0)) { std::cerr << "steps must be > 0\n"; return EXIT_FAILURE; }

    const double H  = static_cast<double>(ny - 2);                        // numer of cells in y-direction
    const double dx = w_phys / H;                                         // lattice spacing (in m)
    const double dt = u_max_lu * dx / u_max_p;                            // time step (in s)
    const double Cu = dx / dt;                                            // conversion factor for velocity
    const double nu_lu = nu_phys * dt / (dx * dx);                        // kinematic viscosity in lattice units
    const double tau = 0.5 + 3.0 * nu_lu;                                 // relaxation time

    const double mu_phys    = rho_phys * nu_phys;                         // dynamic viscosity
    const double dpdx_phys  = 8.0 * mu_phys * u_max_p / (w_phys*w_phys);  // pressure gradient (from analytical Poiseuille solution)
    const double g_phys     = dpdx_phys / rho_phys;                       // body acceleration
    const double force_x = g_phys * (dt * dt / dx);                       // convert to lattice units

    std::cout << "Simulation parameters: nx=" << nx << ", ny=" << ny
              << ", tau=" << tau << ", u_max_lu=" << u_max_lu
              << ", steps=" << steps
              << ", force_x = " << force_x
              << ", dx=" << dx << " [m], dt=" << dt << " [s], Cu=dx/dt=" << Cu << " [m/s per lu]" << std::endl;

    try {
        lbm::Solver solver(nx, ny, tau, force_x, u_max_lu, steps);
        solver.run();
        // Write the resulting velocity profile to disk
        solver.writeVelocityProfile("velocity_profile.csv", dx, Cu);
        std::cout << "Resulting velocity profile successfully written to disk." << std::endl;
    } catch (const std::exception &ex) {
        std::cerr << "Error during simulation: " << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
    return 0;
}
