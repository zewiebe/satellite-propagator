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

Shared physical constants are declared in `include/satellite/constants.h`. The initial values are Earth’s equatorial radius (`6378.137 km`), Earth’s gravitational parameter (`398600.4418 km^3/s^2`), mathematical pi, and constexpr degree/radian conversion factors.

## State vector

`satellite::StateVector` represents the translational state of a spacecraft with a Cartesian position vector in `km` and a Cartesian velocity vector in `km/s`. `satellite::State` remains available as a shorter compatibility name.

For the default Earth two-body model, call `satellite::propagate(initial_state, duration_seconds, step_seconds)`. It uses the Earth gravitational parameter and RK4. Construct `satellite::Propagator` directly when a different gravitational parameter or explicit Euler/RK4 selection is needed.

## Gravity model

`Propagator::acceleration` applies the Earth-centered two-body model

$$\mathbf{a} = -\frac{\mu}{r^3}\mathbf{r}$$

where `mu` is Earth's gravitational parameter and `r` is the position magnitude. The model includes no atmospheric drag, third-body gravity, oblateness, or other perturbations. A position at the central-body origin is rejected because the equation is singular there.

## Euler and RK4 accuracy comparison

The test suite compares both integrators with identical initial conditions: a circular orbit at `7000 km` with velocity `sqrt(mu / r)` in the positive y direction. Each method propagates for one analytical orbital period (`5828.516638 s`) using a `60 s` nominal step; the final step is shortened to end exactly at the period. Position and velocity errors are Euclidean differences from the analytical circular-orbit state at the end of the period. Energy and angular-momentum errors are the absolute changes in specific orbital energy and angular-momentum magnitude from their initial values, expressed as percentages of the initial magnitudes.

| Integrator | Position error (km) | Velocity error (km/s) | Energy error (%) | Angular momentum error (%) |
| --- | ---: | ---: | ---: | ---: |
| Euler | 16055.061167 | 10.999280 | 35.412089 | 23.092463 |
| RK4 | 0.024156 | 0.000026 | 0.000020 | 0.000010 |

For this one-orbit case, RK4 has substantially lower state and invariant errors. These measurements depend on the initial orbit, propagation duration, and step size; they are not a general error bound. The comparison is exercised by `propagator-tests`.

The same test validates the RK4 circular orbit throughout one period using 10-second integration steps and 60-second checkpoints. For the 7000 km orbit, the analytical period is approximately `5828.516638 s`; the test requires the computed period to be within `0.001 s`, sampled radius to stay within `0.001 km`, sampled speed within `0.000001 km/s`, and final position within `0.001 km` of the initial position.

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
