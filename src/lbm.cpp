#include "lbm.h"
#include <cmath>
#include <fstream>
#include <iostream>

namespace lbm {

Solver::Solver(std::size_t nx, std::size_t ny, double tau,
               double force_x, double u_max, std::size_t steps)
    : nx_(nx), ny_(ny), tau_(tau),
      nu_((tau - 0.5) / 3.0),
      force_x_(force_x), u_max_(u_max), steps_(steps),
      f1_(nx * ny * 9, 0.0), f2_(nx * ny * 9, 0.0),
      rho_(nx * ny, 1.0), uX_(nx * ny, 0.0), uY_(nx * ny, 0.0) {
    // Ensure stability requirements.
    if (tau <= 0.5) {
        throw std::invalid_argument("Relaxation time tau must be greater than 0.5");
    }
}

void Solver::initialiseFields() {
    // Initialise density and velocity.
    std::fill(rho_.begin(), rho_.end(), 1.0);
    std::fill(uX_.begin(), uX_.end(), 0.0);
    std::fill(uY_.begin(), uY_.end(), 0.0);
    // Populate the distribution function at equilibrium for each node.
    for (std::size_t y = 0; y < ny_; ++y) {
        for (std::size_t x = 0; x < nx_; ++x) {
            double rho = 1.0;
            double ux = 0.0;
            double uy = 0.0;
            for (std::size_t d = 0; d < 9; ++d) {
                double c_u = dirX_[d] * ux + dirY_[d] * uy;
                double u_sq = ux * ux + uy * uy;
                double feq = weights_[d] * rho * (1.0 + 3.0 * c_u + 4.5 * c_u * c_u - 1.5 * u_sq);
                f1_[idx(x, y, d)] = feq;
            }
        }
    }
}

void Solver::stream(const std::vector<double> &src, std::vector<double> &dest) {
    // Clear destination buffer
    std::fill(dest.begin(), dest.end(), 0.0);
    // Perform streaming for interior points
    for (std::size_t y = 0; y < ny_; ++y) {
        for (std::size_t x = 0; x < nx_; ++x) {
            for (std::size_t d = 0; d < 9; ++d) {
                // Compute source coordinates
                std::size_t prev_x = (x + nx_ - dirX_[d]) % nx_;
                int prev_y_int = static_cast<int>(y) - dirY_[d];
                // Only stream if the source y is inside the domain
                if (prev_y_int >= 0 && prev_y_int < static_cast<int>(ny_)) {
                    std::size_t prev_y = static_cast<std::size_t>(prev_y_int);
                    dest[idx(x, y, d)] = src[idx(prev_x, prev_y, d)];
                }
            }
        }
    }
    // Apply bounce–back on the bottom (y=0) and top (y=ny_-1)
    // Bottom boundary bounce–back
    for (std::size_t x = 0; x < nx_; ++x) {
        std::size_t y = 0;
        for (std::size_t d = 0; d < 9; ++d) {
            if (dirY_[d] == -1) {
                // Find opposite direction
                std::size_t opp = static_cast<std::size_t>(opposite_[d]);
                // Bounce back
                dest[idx(x, y, opp)] = src[idx(x, y, d)];
            }
        }
    }
    // Top boundary bounce–back
    for (std::size_t x = 0; x < nx_; ++x) {
        std::size_t y = ny_ - 1;
        for (std::size_t d = 0; d < 9; ++d) {
            if (dirY_[d] == 1) {
                std::size_t opp = static_cast<std::size_t>(opposite_[d]);
                dest[idx(x, y, opp)] = src[idx(x, y, d)];
            }
        }
    }
}

void Solver::computeMacroscopic(const std::vector<double> &f) {
    // Compute density and velocity from distributions
    for (std::size_t y = 0; y < ny_; ++y) {
        for (std::size_t x = 0; x < nx_; ++x) {
            double rho = 0.0;
            double jx = 0.0;
            double jy = 0.0;
            // Accumulate density and momentum
            for (std::size_t d = 0; d < 9; ++d) {
                double fi = f[idx(x, y, d)];
                rho += fi;
                jx += fi * dirX_[d];
                jy += fi * dirY_[d];
            }
            std::size_t lin = y * nx_ + x;
            rho_[lin] = rho;
            // Compute velocity with Guo forcing correction: momentum/density = velocity
            const double invRho = (rho > 0.0) ? 1.0 / rho : 0.0;
            uX_[lin] = jx * invRho + 0.5 * force_x_;
            uY_[lin] = jy * invRho;
        }
    }
}

void Solver::collide(std::vector<double> &f) {
    // Relax towards equilibrium and apply forcing
    for (std::size_t y = 0; y < ny_; ++y) {
        for (std::size_t x = 0; x < nx_; ++x) {
            std::size_t lin = y * nx_ + x;
            double rho = rho_[lin];
            double ux = uX_[lin];
            double uy = uY_[lin];
            double u_sq = ux * ux + uy * uy;
            for (std::size_t d = 0; d < 9; ++d) {
                // Compute equilibrium distribution
                double ex = static_cast<double>(dirX_[d]);
                double ey = static_cast<double>(dirY_[d]);
                double c_u = ex * ux + ey * uy;
                double feq = weights_[d] * rho * (1.0 + 3.0 * c_u + 4.5 * c_u * c_u - 1.5 * u_sq);

                // Get current population and force parameters
                std::size_t index = idx(x, y, d);
                double fi = f[index];
                double one_minus_half_omega = 1.0 - 1.0 / (2.0 * tau_);
                double Fi = weights_[d] * one_minus_half_omega *
                            ( 3.0 * (ex - ux) + 9.0 * c_u * ex ) * force_x_;

                // Relax towards equilibrium + add force
                f[index] = fi - (fi - feq) / tau_ + Fi;
            }
        }
    }
}

void Solver::run() {
    initialiseFields();
    bool current_is_f1 = true;
    for (std::size_t step = 0; step < steps_; ++step) {
        if (current_is_f1) {
            // Stream from f1_ into f2_
            stream(f1_, f2_);
            // Compute macroscopic variables from f2_
            computeMacroscopic(f2_);
            // Relax distributions in f2_
            collide(f2_);
        } else {
            stream(f2_, f1_);
            computeMacroscopic(f1_);
            collide(f1_);
        }
        current_is_f1 = !current_is_f1;
    }
}

void Solver::writeVelocityProfile(const std::string &filename, double dx, double Cu) const {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Failed to open " << filename << " for writing\n";
        return;
    }
    // Extract the horizontal velocity profile at a chosen x
    std::size_t x = nx_ / 2;
    for (std::size_t y = 0; y < ny_; ++y) {
        double ux = uX_[y * nx_ + x];
        out << (y * dx) << "," << (ux * Cu) << "\n";
    }
    out.close();
}

double Solver::computeError() const {
    // Compute the root mean square error between simulated velocity profile and analytical parabola.
    // Ignore wall points because velocity is zero there.
    double err2 = 0.0;
    std::size_t count = 0;
    std::size_t x = nx_ / 2;
    double H = static_cast<double>(ny_ - 2);
    for (std::size_t y = 1; y + 1 < ny_; ++y) {
        double y_rel = (static_cast<double>(y) - 0.5) / H;
        double u_exact = 4.0 * u_max_ * y_rel * (1.0 - y_rel);
        double u_sim = uX_[y * nx_ + x];
        double diff = u_sim - u_exact;
        err2 += diff * diff;
        ++count;
    }
    return std::sqrt(err2 / static_cast<double>(count));
}

} // namespace lbm
