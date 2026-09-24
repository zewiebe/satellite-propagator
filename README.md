# Satellite Propagator

A standalone C++17 satellite propagation simulator using a fourth-order Runge-Kutta integrator and a two-body Earth model.

## Simulation units

The simulator uses a consistent kilometer-kilogram-second convention:

| Quantity | Unit | Notes |
| --- | --- | --- |
| Distance / position | kilometer (`km`) | Cartesian coordinates use this unit. |
| Velocity | kilometer per second (`km/s`) | State velocity components use this unit. |
| Acceleration | kilometer per second squared (`km/s^2`) | The two-body gravity model returns this unit. |
| Time / duration / step | second (`s`) | Propagation intervals and integration steps use this unit. |
| Angle | radian (`rad`) | Use radians for future attitude or orbital-angle APIs. |
| Mass | kilogram (`kg`) | Use kilograms for future spacecraft or body-mass APIs. |
| Gravitational parameter | cubic kilometer per second squared (`km^3/s^2`) | The constructor's `mu` value is expressed in this unit. |

Angles are measured in radians rather than degrees. Mass is recorded in kilograms, while the current two-body propagator takes the gravitational parameter directly and does not yet expose a mass parameter.

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
