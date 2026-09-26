#pragma once

#include <functional>

namespace satellite {

// Cartesian position vectors are expressed in kilometers (km) when used as
// positions and in kilometers per second (km/s) when used as velocities.
struct Vector3 {
    double x;
    double y;
    double z;

    Vector3 operator+(const Vector3& other) const;
    Vector3 operator-(const Vector3& other) const;
    Vector3 operator*(double scalar) const;
};

// A Cartesian translational state vector containing position and velocity.
// Position is in kilometers (km); velocity is in kilometers per second (km/s).
struct StateVector {
    // Position in kilometers (km).
    Vector3 position;
    // Velocity in kilometers per second (km/s).
    Vector3 velocity;
};

// Compatibility name for callers using the shorter state type.
using State = StateVector;

// Generic first-order state derivative. For translational motion, the state
// derivative is the time rate of change of position and velocity.
struct StateDerivative {
    // Time derivative of the position vector, usually velocity.
    Vector3 position;
    // Time derivative of the velocity vector, usually acceleration.
    Vector3 velocity;

    Vector3 position_derivative() const noexcept { return position; }
    Vector3 velocity_derivative() const noexcept { return velocity; }
};

using StateDerivativeFunction = std::function<StateDerivative(const State&)>;

class Integrator {
public:
    virtual ~Integrator() = default;

    virtual State integrate(const State& initial_state, double duration_seconds,
                            double step_seconds,
                            const StateDerivativeFunction& derivative) const = 0;
};

class EulerIntegrator final : public Integrator {
public:
    State integrate(const State& initial_state, double duration_seconds,
                    double step_seconds,
                    const StateDerivativeFunction& derivative) const override;
};

class RungeKutta4Integrator final : public Integrator {
public:
    State integrate(const State& initial_state, double duration_seconds,
                    double step_seconds,
                    const StateDerivativeFunction& derivative) const override;
};

class Propagator {
public:
    // gravitational_parameter_km3_s2 is mu in km^3/s^2.
    explicit Propagator(double gravitational_parameter_km3_s2);

    // duration_seconds and step_seconds are in seconds (s).
    State propagate(const State& initial_state, double duration_seconds,
                    double step_seconds) const;

    // Simple Euler reference solution for numerical-error comparison.
    State propagate_euler(const State& initial_state, double duration_seconds,
                          double step_seconds) const;

    // Fourth-order Runge-Kutta integration for the production propagator.
    State propagate_rk4(const State& initial_state, double duration_seconds,
                        double step_seconds) const;

    // Generic derivative-based propagation using an externally supplied state
    // rate model. This keeps the numerical integrator separate from the physics.
    template <typename DerivativeFunction>
    State propagate(const State& initial_state, double duration_seconds,
                    double step_seconds,
                    const DerivativeFunction& derivative) const {
        return RungeKutta4Integrator{}.integrate(initial_state, duration_seconds,
                                                step_seconds,
                                                StateDerivativeFunction{derivative});
    }

    // Returns Earth-centered two-body acceleration in km/s^2. The supplied
    // gravitational parameter is the only force parameter; no perturbations
    // are modeled.
    Vector3 acceleration(const Vector3& position) const;

    // State time derivative for the two-body gravity model.
    StateDerivative derivative(const State& state) const;

private:
    double gravitational_parameter_;
};

// Propagates with the Earth gravitational parameter and the default RK4
// integrator. Duration and timestep are in seconds.
StateVector propagate(const StateVector& initial_state, double duration_seconds,
                      double step_seconds);

}
