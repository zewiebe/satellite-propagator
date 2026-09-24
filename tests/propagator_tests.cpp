#include "satellite/propagator.h"
#include "satellite/constants.h"
#include "math/vector3.h"

#include <cmath>
#include <iostream>

namespace {
bool close(double actual, double expected, double tolerance) {
    return std::abs(actual - expected) <= tolerance;
}

bool close_vector(const math::Vector3& actual, const math::Vector3& expected,
                  double tolerance) {
    return close(actual.x, expected.x, tolerance) &&
           close(actual.y, expected.y, tolerance) &&
           close(actual.z, expected.z, tolerance);
}

bool close_vector(const satellite::Vector3& actual,
            const satellite::Vector3& expected, double tolerance) {
    return close(actual.x, expected.x, tolerance) &&
        close(actual.y, expected.y, tolerance) &&
        close(actual.z, expected.z, tolerance);
}

double magnitude(const satellite::Vector3& vector) {
    return std::sqrt(vector.x * vector.x + vector.y * vector.y +
                     vector.z * vector.z);
}
}

int main() {
    constexpr double tolerance = 1e-12;

    if (!close(satellite::constants::earth_radius_km, 6378.137, tolerance) ||
        !close(satellite::constants::earth_gravitational_parameter_km3_s2,
               398600.4418, tolerance) ||
        !close(satellite::constants::pi, 3.141592653589793, tolerance) ||
        !close(satellite::constants::degrees_to_radians * 180.0,
               satellite::constants::pi, tolerance) ||
        !close(satellite::constants::radians_to_degrees *
                   satellite::constants::pi,
               180.0, tolerance)) {
        std::cerr << "Unexpected physical constants\n";
        return 1;
    }

    const math::Vector3 first{1.0, -2.0, 3.0};
    const math::Vector3 second{-4.0, 5.0, 6.0};
    if (!close_vector(first + second, {-3.0, 3.0, 9.0}, tolerance)) {
        std::cerr << "Unexpected vector addition\n";
        return 1;
    }
    if (!close_vector(first - second, {5.0, -7.0, -3.0}, tolerance)) {
        std::cerr << "Unexpected vector subtraction\n";
        return 1;
    }

    const math::Vector3 magnitude_vector{3.0, 4.0, 12.0};
    if (!close(magnitude_vector.magnitude(), 13.0, tolerance) ||
        !close(magnitude_vector.magnitude_squared(), 169.0, tolerance)) {
        std::cerr << "Unexpected vector magnitude\n";
        return 1;
    }

    const math::Vector3 normal = math::Vector3{3.0, 4.0, 0.0}.normalized();
    if (!close_vector(normal, {0.6, 0.8, 0.0}, tolerance) ||
        !close(normal.magnitude(), 1.0, tolerance)) {
        std::cerr << "Unexpected vector normalization\n";
        return 1;
    }

    const math::Vector3 zero{};
    if (!close_vector(zero.normalized(), zero, tolerance)) {
        std::cerr << "Zero vector normalization should remain zero\n";
        return 1;
    }

    if (!close(first.dot(second), 4.0, tolerance) ||
        !close(math::Vector3{1.0, 0.0, 0.0}.dot({0.0, 1.0, 0.0}), 0.0,
               tolerance)) {
        std::cerr << "Unexpected vector dot product\n";
        return 1;
    }

    if (!close_vector(math::Vector3{1.0, 0.0, 0.0}.cross({0.0, 1.0, 0.0}),
                      {0.0, 0.0, 1.0}, tolerance) ||
        !close_vector(first.cross(first), zero, tolerance) ||
        !close_vector(zero.cross(first), zero, tolerance)) {
        std::cerr << "Unexpected vector cross product\n";
        return 1;
    }

    satellite::Propagator propagator(
        satellite::constants::earth_gravitational_parameter_km3_s2);
    const satellite::Vector3 acceleration = propagator.acceleration({7000.0, 0.0, 0.0});

    if (!close(acceleration.x,
               -satellite::constants::earth_gravitational_parameter_km3_s2 /
                   (7000.0 * 7000.0),
               1e-12) ||
        !close(acceleration.y, 0.0, 1e-12) ||
        !close(acceleration.z, 0.0, 1e-12)) {
        std::cerr << "Unexpected gravitational acceleration\n";
        return 1;
    }

    const satellite::Vector3 three_dimensional_position{3000.0, 4000.0, 12000.0};
    const satellite::Vector3 three_dimensional_acceleration =
        propagator.acceleration(three_dimensional_position);
    const double radius = std::sqrt(
        three_dimensional_position.x * three_dimensional_position.x +
        three_dimensional_position.y * three_dimensional_position.y +
        three_dimensional_position.z * three_dimensional_position.z);
    const double expected_scale =
        -satellite::constants::earth_gravitational_parameter_km3_s2 /
        (radius * radius * radius);
    if (!close(three_dimensional_acceleration.x,
               three_dimensional_position.x * expected_scale, tolerance) ||
        !close(three_dimensional_acceleration.y,
               three_dimensional_position.y * expected_scale, tolerance) ||
        !close(three_dimensional_acceleration.z,
               three_dimensional_position.z * expected_scale, tolerance)) {
        std::cerr << "Gravity did not follow the two-body position vector\n";
        return 1;
    }

    const auto check_gravity_at_radius = [&](const satellite::Vector3& position,
                                             double expected_radius,
                                             const char* case_name) {
        const satellite::Vector3 acceleration = propagator.acceleration(position);
        const double acceleration_magnitude = magnitude(acceleration);
        const double expected_magnitude =
            satellite::constants::earth_gravitational_parameter_km3_s2 /
            (expected_radius * expected_radius);
        const double expected_scale = -
            satellite::constants::earth_gravitational_parameter_km3_s2 /
            (expected_radius * expected_radius * expected_radius);

        if (!close(magnitude(position), expected_radius, tolerance) ||
            !close(acceleration_magnitude, expected_magnitude, tolerance) ||
            !close_vector(acceleration, position * expected_scale, tolerance)) {
            std::cerr << "Unexpected gravitational acceleration at "
                      << case_name << "\n";
            return false;
        }
        return true;
    };

    if (!check_gravity_at_radius(
            {satellite::constants::earth_radius_km, 0.0, 0.0},
            satellite::constants::earth_radius_km, "Earth surface") ||
        !check_gravity_at_radius(
            {satellite::constants::earth_radius_km + 400.0, 0.0, 0.0},
            satellite::constants::earth_radius_km + 400.0, "400 km altitude") ||
        !check_gravity_at_radius({0.0, 0.0, 10000.0}, 10000.0,
                                  "higher altitude")) {
        return 1;
    }

    const satellite::StateVector initial_state{{7000.0, 0.0, 0.0},
                                               {0.0, 7.54605329, 0.0}};
    const satellite::StateVector unchanged_state =
        propagator.propagate(initial_state, 0.0, 1.0);
    if (!close(unchanged_state.position.x, initial_state.position.x, tolerance) ||
        !close(unchanged_state.position.y, initial_state.position.y, tolerance) ||
        !close(unchanged_state.position.z, initial_state.position.z, tolerance) ||
        !close(unchanged_state.velocity.x, initial_state.velocity.x, tolerance) ||
        !close(unchanged_state.velocity.y, initial_state.velocity.y, tolerance) ||
        !close(unchanged_state.velocity.z, initial_state.velocity.z, tolerance)) {
        std::cerr << "State vector components were not preserved\n";
        return 1;
    }

    const satellite::State after_one_minute = propagator.propagate(initial_state, 60.0, 1.0);
    if (!(after_one_minute.position.x < initial_state.position.x &&
          after_one_minute.position.y > initial_state.position.y)) {
        std::cerr << "Orbit did not advance in the expected direction\n";
        return 1;
    }

    return 0;
}
