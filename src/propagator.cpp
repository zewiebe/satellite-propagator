#include "satellite/propagator.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace satellite {

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

    const double radius_cubed = radius_squared * std::sqrt(radius_squared);
    return position * (-gravitational_parameter_ / radius_cubed);
}

State Propagator::propagate(const State& initial_state, double duration_seconds,
                            double step_seconds) const {
    if (duration_seconds < 0.0 || step_seconds <= 0.0) {
        throw std::invalid_argument("duration must be non-negative and step must be positive");
    }

    State state = initial_state;
    double elapsed = 0.0;
    while (elapsed < duration_seconds) {
        const double step = std::min(step_seconds, duration_seconds - elapsed);
        const Vector3 k1_position = state.velocity;
        const Vector3 k1_velocity = acceleration(state.position);

        const State midpoint_one{state.position + k1_position * (step * 0.5),
                                 state.velocity + k1_velocity * (step * 0.5)};
        const Vector3 k2_position = midpoint_one.velocity;
        const Vector3 k2_velocity = acceleration(midpoint_one.position);

        const State midpoint_two{state.position + k2_position * (step * 0.5),
                                 state.velocity + k2_velocity * (step * 0.5)};
        const Vector3 k3_position = midpoint_two.velocity;
        const Vector3 k3_velocity = acceleration(midpoint_two.position);

        const State endpoint{state.position + k3_position * step,
                             state.velocity + k3_velocity * step};
        const Vector3 k4_position = endpoint.velocity;
        const Vector3 k4_velocity = acceleration(endpoint.position);

        state.position = state.position +
                         (k1_position + k2_position * 2.0 + k3_position * 2.0 + k4_position) * (step / 6.0);
        state.velocity = state.velocity +
                         (k1_velocity + k2_velocity * 2.0 + k3_velocity * 2.0 + k4_velocity) * (step / 6.0);
        elapsed += step;
    }

    return state;
}

}
