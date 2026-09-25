#include "satellite/propagator.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace satellite {

namespace {

void validate_duration_and_step(double duration_seconds, double step_seconds) {
    if (duration_seconds < 0.0 || step_seconds <= 0.0) {
        throw std::invalid_argument(
            "duration must be non-negative and step must be positive");
    }
}

State integrate_with_step(const State& initial_state, double duration_seconds,
                          double step_seconds,
                          const StateDerivativeFunction& derivative,
                          bool use_euler) {
    if (!derivative) {
        throw std::invalid_argument("state derivative must be supplied");
    }

    validate_duration_and_step(duration_seconds, step_seconds);

    State state = initial_state;
    double elapsed = 0.0;
    while (elapsed < duration_seconds) {
        const double step = std::min(step_seconds, duration_seconds - elapsed);

        if (use_euler) {
            const StateDerivative derivative_value = derivative(state);
            state.position = state.position + derivative_value.position * step;
            state.velocity = state.velocity + derivative_value.velocity * step;
        } else {
            const StateDerivative k1 = derivative(state);
            const State midpoint_one{state.position + k1.position * (step * 0.5),
                                     state.velocity + k1.velocity * (step * 0.5)};
            const StateDerivative k2 = derivative(midpoint_one);
            const State midpoint_two{state.position + k2.position * (step * 0.5),
                                     state.velocity + k2.velocity * (step * 0.5)};
            const StateDerivative k3 = derivative(midpoint_two);
            const State endpoint{state.position + k3.position * step,
                                 state.velocity + k3.velocity * step};
            const StateDerivative k4 = derivative(endpoint);

            state.position = state.position +
                             (k1.position + k2.position * 2.0 + k3.position * 2.0 +
                              k4.position) * (step / 6.0);
            state.velocity = state.velocity +
                             (k1.velocity + k2.velocity * 2.0 + k3.velocity * 2.0 +
                              k4.velocity) * (step / 6.0);
        }

        elapsed += step;
    }

    return state;
}

} // namespace

Vector3 Vector3::operator+(const Vector3& other) const {
    return {x + other.x, y + other.y, z + other.z};
}

Vector3 Vector3::operator-(const Vector3& other) const {
    return {x - other.x, y - other.y, z - other.z};
}

Vector3 Vector3::operator*(double scalar) const {
    return {x * scalar, y * scalar, z * scalar};
}

Propagator::Propagator(double gravitational_parameter_km3_s2)
    : gravitational_parameter_(gravitational_parameter_km3_s2) {
    if (gravitational_parameter_ <= 0.0) {
        throw std::invalid_argument("gravitational parameter must be positive");
    }
}

Vector3 Propagator::acceleration(const Vector3& position) const {
    const double radius_squared = position.x * position.x +
                                  position.y * position.y +
                                  position.z * position.z;
    if (radius_squared == 0.0) {
        throw std::invalid_argument("position cannot be at the central body origin");
    }

    const double radius = std::sqrt(radius_squared);
    const double radius_cubed = radius_squared * radius;
    return position * (-gravitational_parameter_ / radius_cubed);
}

StateDerivative Propagator::derivative(const State& state) const {
    return {state.velocity, acceleration(state.position)};
}

State Propagator::propagate(const State& initial_state, double duration_seconds,
                            double step_seconds) const {
    return propagate_rk4(initial_state, duration_seconds, step_seconds);
}

State Propagator::propagate_euler(const State& initial_state,
                                 double duration_seconds,
                                 double step_seconds) const {
    return EulerIntegrator{}.integrate(initial_state, duration_seconds,
                                       step_seconds,
                                       [this](const State& state) {
                                           return derivative(state);
                                       });
}

State Propagator::propagate_rk4(const State& initial_state,
                               double duration_seconds,
                               double step_seconds) const {
    return RungeKutta4Integrator{}.integrate(initial_state, duration_seconds,
                                            step_seconds,
                                            [this](const State& state) {
                                                return derivative(state);
                                            });
}

State EulerIntegrator::integrate(const State& initial_state,
                                double duration_seconds, double step_seconds,
                                const StateDerivativeFunction& derivative) const {
    return integrate_with_step(initial_state, duration_seconds, step_seconds,
                              derivative, true);
}

State RungeKutta4Integrator::integrate(const State& initial_state,
                                      double duration_seconds,
                                      double step_seconds,
                                      const StateDerivativeFunction& derivative) const {
    return integrate_with_step(initial_state, duration_seconds, step_seconds,
                              derivative, false);
}

}
