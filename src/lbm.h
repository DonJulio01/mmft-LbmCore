#pragma once

#include <vector>
#include <cstddef>
#include <string>

namespace lbm {

/** D2Q9 LBM (BGK) for 2D Poiseuille. No-slip at y=0,ny-1; periodic in x; body force in x. */
class Solver {
public:
    /**
     * Construct a lattice Boltzmann solver.
     *
     * @param nx Number of lattice nodes in the x (streamwise) direction.
     * @param ny Number of lattice nodes in the y (wall–normal) direction.
     * @param tau Relaxation time controlling the kinematic viscosity via
     *            nu = (tau-0.5)/3.  Values must satisfy tau > 0.5 for
     *            stability.
     * @param force_x Constant body force in lattice units acting along x.
     * @param u_max Desired maximum velocity used for analytical reference.
     * @param steps Number of time steps to simulate.
     */
    Solver(std::size_t nx, std::size_t ny, double tau,
           double force_x, double u_max, std::size_t steps);

    /**
     * Run the time stepping loop.  Distributions are initialised to
     * equilibrium at rest before integration.
     */
    void run();

    /** Write the velocity profile to csv. */
    void writeVelocityProfile(const std::string &filename) const;

    /**
     * Compute an L2–norm error between the simulated velocity profile and
     * the analytical Poiseuille profile.  The analytical profile is
     * defined as u(y) = 4 * u_max * (y/H) * (1 - y/H), where H is the
     * number of lattice spacings across the channel minus one.  Only the
     * x–component of the velocity is compared.
     *
     * @return The root mean square error of u_x compared to the analytic
     *         solution.
     */
    double computeError() const;

private:
    std::size_t nx_;
    std::size_t ny_;
    double tau_;
    double nu_;
    double force_x_;
    double u_max_;
    std::size_t steps_;

    // Discrete velocity directions for the D2Q9 scheme
    inline static constexpr int dirX_[9]  = {0, 1, 0, -1,  0, 1, -1, -1,  1};
    inline static constexpr int dirY_[9]  = {0, 0, 1,  0, -1, 1,  1, -1, -1};
    inline static constexpr int opposite_[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    inline static constexpr double weights_[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

    // Distribution functions (two buffers)
    std::vector<double> f1_;
    std::vector<double> f2_;
    // Macroscopic fields.
    std::vector<double> rho_;
    std::vector<double> uX_;
    std::vector<double> uY_;

    // Compute a linear index into the flattened f1_ / f2_ vectors.
    inline std::size_t idx(std::size_t x, std::size_t y, std::size_t d) const {
        return (y * nx_ + x) * 9 + d;
    }

    /**
     * Initialise the simulation state.  Density is set to unity and
     * velocities are zero so that the equilibrium distribution can be
     * computed accordingly.
     */
    void initialiseFields();

    /**
     * Perform the streaming step.  The source and destination buffers
     * must be distinct.  Periodicity along x is handled implicitly and
     * bounce–back is applied for nodes on the top and bottom boundaries.
     *
     * @param src Source distribution array to stream from.
     * @param dest Destination distribution array to stream into.
     */
    void stream(const std::vector<double> &src, std::vector<double> &dest);

    /**
     * Compute density and velocity fields from the distribution function.
     *
     * @param f The distribution array to sample.
     */
    void computeMacroscopic(const std::vector<double> &f);

    /**
     * Perform the collision step using the BGK operator and include a
     * forcing term in the x direction.  The distribution is updated in
     * place.
     *
     * @param f Distribution array to relax towards equilibrium and apply forcing.
     */
    void collide(std::vector<double> &f);
};

} // namespace lbm
