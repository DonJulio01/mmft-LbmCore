#pragma once

#include <vector>
#include <cstddef>
#include <string>

namespace lbm {

/** D2Q9 BGK solver for 2D Poiseuille. No-slip at y=0,ny-1. Zou/He inlet (u_x) and outlet (rho). */
class Solver {
public:
    /**
     * Create LBM Solver.
     * @param nx Lattice size in x.
     * @param ny Lattice size in y.
     * @param tau Relaxation time (tau>0.5).
     * @param u_max Peak parabolic inlet speed.
     * @param steps Time steps to run.
     * @param rho_out_lu Outlet density in lattice units.
     */
    Solver(std::size_t nx, std::size_t ny, double tau, double u_max, std::size_t steps, double rho_out_lu = 1.0);

    /**
     * Run the time stepping loop.
     */
    void run();

    /** Write u_x(y) at x=nx/2 to CSV. dx is spacing; Cu converts to physical if >0. */
    void writeVelocityProfile(const std::string &filename, double dx, double Cu) const;

    /**
     * L2 error vs u(y)=4*u_max*(y/H)*(1-y/H) using u_x.
     * @return RMS error.
     */
    double computeError() const;

private:
    std::size_t nx_;
    std::size_t ny_;
    double tau_;
    double nu_;
    double u_max_;
    std::size_t steps_;
    double rho_out_;

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
     * velocities are zero.
     */
    void initialiseFields();

    /**
     * Pull streaming with bounce-back; Zou/He inlet (u_x) and outlet (rho).
     * @param src Source distributions.
     * @param dest Destination distributions.
     */
    void stream(const std::vector<double> &src, std::vector<double> &dest);

    /**
     * Compute density and velocity fields from the distribution function.
     *
     * @param f The distribution array to sample.
     */
    void computeMacroscopic(const std::vector<double> &f);

    /**
     * BGK collision in place.
     * @param f Distributions to relax.
     */
    void collide(std::vector<double> &f);
};

} // namespace lbm
