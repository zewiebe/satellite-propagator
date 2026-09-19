# Satellite Propagator

A standalone C++17 satellite propagation simulator using a fourth-order Runge-Kutta integrator and a two-body Earth model. Positions are expressed in kilometers and velocities in kilometers per second.

## Build

Requirements: CMake 3.16 or newer and a C++17 compiler.

```powershell
cmake --preset tdm-gcc-mingw
cmake --build --preset tdm-gcc-mingw-release
ctest --preset tdm-gcc-mingw-release
```

The old `build` directory contains an NMake cache. Do not reuse it with MinGW; use the `build-mingw` preset directory above, or delete `build/CMakeCache.txt` and `build/CMakeFiles` before configuring `build` again.

## Run

The executable propagates a 7000 km circular-orbit starting state for one hour by default and prints the final state:

```powershell
.\build-mingw\satellite-propagator.exe
```

Pass a duration and integration step in seconds to override the defaults:

```powershell
.\build-mingw\satellite-propagator.exe 600 1
```

The propagation library is available through `include/satellite/propagator.h` for use by other C++ programs.
