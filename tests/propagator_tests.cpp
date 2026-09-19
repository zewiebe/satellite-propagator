#include "satellite/propagator.h"

#include <cmath>
#include <iostream>

namespace {
constexpr double earth_mu_km3_s2 = 398600.4418;

bool close(double actual, double expected, double tolerance) {
    return std::abs(actual - expected) <= tolerance;
}
}

int main() {
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
