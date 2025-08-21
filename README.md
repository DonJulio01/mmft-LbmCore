# mmft-LbmCore
A lightweight 2D LBM (D2Q9, BGK) solver focused on plane Poiseuille flow.

BCs: no-slip walls (bounce-back) at y=0 and y=ny-1, inlet/outlet boundary conditions in x.     
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
cmake -B build -A x64 -G "Visual Studio 17 2022"
cmake --build build --config Release
```

The executable is build/Release/MMFTLBM.exe

-------------------------------------------------------------------------------
3) Run
-------------------------------------------------------------------------------
**CLI:** MMFTLBM [w_phys] [nu_phys] [rho_phys] [u_max_phys] [u_max_lu] [nx] [ny] [steps] 

Arguments:

  w_phys     — channel width [m]  
  nu_phys    — kinematic viscosity [m²/s]  
  rho_phys   — density [kg/m³]  
  u_max_phys — target peak velocity [m/s]  
  u_max_lu   — target peak velocity in lattice units (for numerical stability)  
  nx, ny     — lattice size in x (streamwise) and y (wall-normal) directions  
  steps      — number of time steps

Example:

```.\build\Release\MMFTLBM.exe 0.01 1.0e-5 1000 0.05 0.05 100 50 20000```

Output:
  velocity_profile.csv  
  Two columns: y [m] (vertical position), u_x(y) [m/s] at x = nx/2 (streamwise flow speed)

Analytical profile:
  u_exact(y) = 4 * u_max_phys * (y/w) * (1 – y/w),  y ∈ [0, w]

-------------------------------------------------------------------------------
4) Tests
-------------------------------------------------------------------------------
Build & run tests

Windows:
```
cmake -B build -A x64 -G "Visual Studio 17 2022" -DTEST=ON
cmake --build build --config Release
ctest --test-dir build -C Release --verbose
```

The test constructs a Poiseuille case and asserts the RMS error between numerical
and analytical profiles is below a tolerance (e.g., 1e-3).
