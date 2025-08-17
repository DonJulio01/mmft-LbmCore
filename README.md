# mmft-LbmCore
A lightweight 2D LBM (D2Q9, BGK) solver focused on plane Poiseuille flow.

BCs: no-slip walls (bounce-back) at y=0 and y=ny-1, periodic in x
Driving: constant body force F_x (Guo forcing) to emulate a pressure gradient
Validation: compares simulated velocity profile with analytical parabola

-------------------------------------------------------------------------------
1) Requirements
-------------------------------------------------------------------------------
- C++17 compiler (MSVC 2022 / clang / gcc)
- CMake ≥ 3.20
- (optional) GoogleTest

-------------------------------------------------------------------------------
2) Build
-------------------------------------------------------------------------------
Windows
```
cd mmft-LbmCore
cmake -B build -A x64 -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

The executable is build/Release/MMFTLBM.exe

-------------------------------------------------------------------------------
3) Run
-------------------------------------------------------------------------------
CLI:
  MMFTLBM [nx] [ny] [tau] [u_max] [steps]

Defaults:
  nx=100  ny=50  tau=0.7  u_max=0.05  steps=10000

Example:
``.\build\Release\MMFTLBM.exe 100 50 0.7 0.05 20000``

Output
  velocity_profile.csv — two columns: y_index, u_x(y) at x = nx/2

-------------------------------------------------------------------------------
4) Parameters
-------------------------------------------------------------------------------
- nx, ny — lattice size in x (streamwise) and y (wall-normal) direction.
- tau — BGK relaxation time. Kinematic viscosity in lattice units: nu = (tau - 0.5)/3.
        Must satisfy tau > 0.5 for stability.
- u_max — target peak velocity of the analytical Poiseuille parabola.
- steps — number of time steps.

Body force computed internally:
  H   = ny - 1
  nu  = (tau - 0.5)/3
  F_x = 8 * nu * u_max / H^2

Analytical profile:
  u_exact(y) = 4 * u_max * (y/H) * (1 - y/H),  y = 0..H

-------------------------------------------------------------------------------
5) Tests
-------------------------------------------------------------------------------
Build & run tests

Windows:
```
cmake -B build -G "Visual Studio 17 2022" -A x64 -DTEST=ON
cmake --build build --config Release
ctest --test-dir build -C Release --verbose
```

The test constructs a Poiseuille case and asserts the RMS error between numerical
and analytical profiles is below a tolerance (e.g., 1e-3).
