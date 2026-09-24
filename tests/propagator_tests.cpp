#include "satellite/propagator.h"
#include "math/vector3.h"

#include <cmath>
#include <iostream>

namespace {
constexpr double earth_mu_km3_s2 = 398600.4418;

bool close(double actual, double expected, double tolerance) {
    return std::abs(actual - expected) <= tolerance;
}

bool close_vector(const math::Vector3& actual, const math::Vector3& expected,
                  double tolerance) {
    return close(actual.x, expected.x, tolerance) &&
           close(actual.y, expected.y, tolerance) &&
           close(actual.z, expected.z, tolerance);
}
}

int main() {
    constexpr double tolerance = 1e-12;

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

    satellite::Propagator propagator(earth_mu_km3_s2);
    const satellite::Vector3 acceleration = propagator.acceleration({7000.0, 0.0, 0.0});

    if (!close(acceleration.x, -earth_mu_km3_s2 / (7000.0 * 7000.0), 1e-12) ||
        !close(acceleration.y, 0.0, 1e-12) ||
        !close(acceleration.z, 0.0, 1e-12)) {
        std::cerr << "Unexpected gravitational acceleration\n";
        return 1;
    }

    const satellite::State initial_state{{7000.0, 0.0, 0.0}, {0.0, 7.54605329, 0.0}};
    const satellite::State after_one_minute = propagator.propagate(initial_state, 60.0, 1.0);
    if (!(after_one_minute.position.x < initial_state.position.x &&
          after_one_minute.position.y > initial_state.position.y)) {
        std::cerr << "Orbit did not advance in the expected direction\n";
        return 1;
    }

    return 0;
}
