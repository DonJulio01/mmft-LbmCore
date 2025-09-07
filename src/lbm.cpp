#include "lbm.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>

namespace lbm {

Solver::Solver(std::size_t nx, std::size_t ny, double tau, double u_max_lu, std::size_t steps, double rho_out_lu):
      nx_(nx), ny_(ny), tau_(tau),
      nu_((tau - 0.5) / 3.0),
      u_max_(u_max_lu), steps_(steps), rho_out_(rho_out_lu),
      f1_(nx * ny * 9, 0.0), f2_(nx * ny * 9, 0.0),
      rho_(nx * ny, 1.0), uX_(nx * ny, 0.0), uY_(nx * ny, 0.0) {
    // Ensure stability requirements.
    if (tau <= 0.5 || tau >= 2.0) {
        throw std::invalid_argument("Relaxation time tau must be greater than 0.5");
    }
}

void Solver::initialiseFields() {
    // Initialise density and velocity.
    std::fill(rho_.begin(), rho_.end(), rho_out_);
    std::fill(uX_.begin(), uX_.end(), 0.0);
    std::fill(uY_.begin(), uY_.end(), 0.0);
    // Populate the distribution function at equilibrium for each node.
    for (std::size_t y = 0; y < ny_; ++y) {
        for (std::size_t x = 0; x < nx_; ++x) {
            double rho = rho_out_;
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
    f2_ = f1_;
}

void Solver::stream(const std::vector<double> &src, std::vector<double> &dest) {
    // Clear destination buffer
    std::fill(dest.begin(), dest.end(), 0.0);
	// Map directions
	int E=1,N=2,W=3,S=4,NE=5,NW=6,SW=7,SE=8;

    // Perform streaming for interior points
    for (std::size_t y = 0; y < ny_; ++y) {
        for (std::size_t x = 0; x < nx_; ++x) {
            for (std::size_t d = 0; d < 9; ++d) {
                // Compute source coordinates
                int prev_x_int = static_cast<int>(x) - dirX_[d];
                int prev_y_int = static_cast<int>(y) - dirY_[d];
                // Only stream if the source is inside the domain
                if (prev_x_int >= 0 && prev_x_int < static_cast<int>(nx_) && prev_y_int >= 0 && prev_y_int < static_cast<int>(ny_)) {
                    std::size_t prev_x = static_cast<std::size_t>(prev_x_int);
                    std::size_t prev_y = static_cast<std::size_t>(prev_y_int);
                    dest[idx(x, y, d)] = src[idx(prev_x, prev_y, d)];
                }
            }
        }
    }

    // No-slip bounce-back
    std::size_t yb = 1, yt = ny_ - 2;
    for (std::size_t x = 0; x < nx_; ++x) {
        for (int d = 0; d < 9; ++d) if (dirY_[d] > 0) dest[idx(x, yb, d)] = src[idx(x, yb, opposite_[d])]; // Bottom wall
        for (int d = 0; d < 9; ++d) if (dirY_[d] < 0) dest[idx(x, yt, d)] = src[idx(x, yt, opposite_[d])]; // Top wall
    }

    // [See 5.3.4.4 in The Lattice Boltzmann Method: Principles and Practice (Timm Krüger et al., 2017)]
    // Inlet (x=0): impose velocity profile (Zou/He), u_y=0, u_x = parabola with peak u_max_
    if (nx_ >= 2 && ny_ >= 3) {
        double H = static_cast<double>(ny_ - 2);
        for (std::size_t y = 1; y + 1 < ny_; ++y) {
            double y_rel = (static_cast<double>(y) - 0.5) / H;
            double ux = 4.0 * u_max_ * y_rel * (1.0 - y_rel);
            // Known: f0,f2,f3,f4,f6,f7 at (0,y) after streaming/bounce-back
            std::size_t x0 = 0;
            double f0 = dest[idx(x0, y, 0)];
            double fN = dest[idx(x0, y, N)];
            double fS = dest[idx(x0, y, S)];
            double fW = dest[idx(x0, y, W)];
            double fNW = dest[idx(x0, y, NW)];
            double fSW = dest[idx(x0, y, SW)];
            double rho = (f0 + fN + fS + 2.0 * (fW + fNW + fSW)) / (1.0 - ux);
            dest[idx(x0, y, E)] = fW + (2.0 / 3.0) * rho * ux;
            dest[idx(x0, y, NE)] = fSW + 0.5 * (fS - fN) + (1.0 / 6.0) * rho * ux;
            dest[idx(x0, y, SE)] = fNW + 0.5 * (fN - fS) + (1.0 / 6.0) * rho * ux;
        }
    }
    // Outlet (x=nx_-1): impose pressure (density) with Zou/He, use rho_out_, u_y=0
    if (nx_ >= 2 && ny_ >= 3) {
        const double rho_out = rho_out_;
        std::size_t xe = nx_ - 1;
        for (std::size_t y = 1; y + 1 < ny_; ++y) {
            // Known: f0,f1,f2,f4,f5,f8 at (xe,y)
            double f0 = dest[idx(xe, y, 0)];
            double fN = dest[idx(xe, y, N)];
            double fS = dest[idx(xe, y, S)];
            double fE = dest[idx(xe, y, E)];
            double fNE = dest[idx(xe, y, NE)];
            double fSE = dest[idx(xe, y, SE)];
            double ux = (f0 + fN + fS + 2.0 * (fE + fNE + fSE)) / rho_out - 1.0;
            dest[idx(xe, y, W)] = fE - (2.0 / 3.0) * rho_out * ux;
            dest[idx(xe, y, NW)] = fSE - 0.5 * (fN - fS) - (1.0 / 6.0) * rho_out * ux;
            dest[idx(xe, y, SW)] = fNE - 0.5 * (fS - fN) - (1.0 / 6.0) * rho_out * ux;
        }
    }
}

void Solver::computeMacroscopic(const std::vector<double> &f) {
    // Compute density and velocity from distributions
    for (std::size_t y = 1; y + 1 < ny_; ++y) {
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
            // Compute velocity: momentum/density = velocity
            double invRho = 1.0 / rho;
            uX_[lin] = jx * invRho;
            uY_[lin] = jy * invRho;
        }
    }
}

void Solver::collide(std::vector<double> &f) {
    // Relax towards equilibrium
    for (std::size_t y = 1; y + 1 < ny_; ++y) {
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

                // Relax towards equilibrium
                f[index] = fi - (fi - feq) / tau_;
            }
        }
    }
}

void Solver::run() {
    initialiseFields();
    std::vector<double>* A = &f1_;
    std::vector<double>* B = &f2_;
    for (std::size_t step = 0; step < steps_; ++step) {
        computeMacroscopic(*A);         // from current populations
        collide(*A);                 // relax in-place
        stream(*A, *B);     // pull streaming + bounce-back + inlet/outlet
        std::swap(A, B);
    }
    computeMacroscopic(*A);
}

void Solver::writeVelocityProfile(const std::string &filename, double dx, double Cu) const {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Failed to open " << filename << " for writing\n";
        return;
    }
    // Extract the horizontal velocity profile at a chosen x
    std::size_t x = nx_ / 2;
    bool phys = std::isfinite(Cu) && Cu > 0.0;
    for (std::size_t y = 1; y + 1 < ny_; ++y) {
        double ux = uX_[y * nx_ + x];
        if (!std::isfinite(ux)) continue;
        double y_phys = (static_cast<double>(y) - 0.5) * dx;
        double v = phys ? ux * Cu : ux;
        out << y_phys << "," << v << "\n";
    }
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
