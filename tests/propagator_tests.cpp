#include "satellite/propagator.h"
#include "satellite/constants.h"
#include "math/vector3.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
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

satellite::Vector3 cross(const satellite::Vector3& first,
                         const satellite::Vector3& second) {
    return {first.y * second.z - first.z * second.y,
            first.z * second.x - first.x * second.z,
            first.x * second.y - first.y * second.x};
}

double specific_energy(const satellite::State& state, double gravitational_parameter) {
    const double speed_squared = state.velocity.x * state.velocity.x +
                                 state.velocity.y * state.velocity.y +
                                 state.velocity.z * state.velocity.z;
    return 0.5 * speed_squared -
           gravitational_parameter / magnitude(state.position);
}

double angular_momentum_magnitude(const satellite::State& state) {
    return magnitude(cross(state.position, state.velocity));
}

struct AccuracyErrors {
    double position_error_km;
    double velocity_error_km_s;
    double energy_error_percent;
    double angular_momentum_error_percent;
};
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

    const satellite::StateVector default_propagated_state =
        satellite::propagate(initial_state, 60.0, 1.0);
    if (!close_vector(default_propagated_state.position,
                      after_one_minute.position, tolerance) ||
        !close_vector(default_propagated_state.velocity,
                      after_one_minute.velocity, tolerance)) {
        std::cerr << "Default propagation interface did not use Earth RK4 propagation\n";
        return 1;
    }

    const auto constant_velocity_derivative =
        [](const satellite::State& state) {
            return satellite::StateDerivative{state.velocity, {0.0, 0.0, 0.0}};
        };
    const satellite::State euler_state =
        satellite::EulerIntegrator{}.integrate({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
                                              1.0, 0.5,
                                              constant_velocity_derivative);
    if (!close(euler_state.position.x, 1.0, tolerance) ||
        !close(euler_state.position.y, 0.0, tolerance) ||
        !close(euler_state.position.z, 0.0, tolerance) ||
        !close(euler_state.velocity.x, 1.0, tolerance) ||
        !close(euler_state.velocity.y, 0.0, tolerance) ||
        !close(euler_state.velocity.z, 0.0, tolerance)) {
        std::cerr << "Euler integrator did not respect the generic state derivative\n";
        return 1;
    }

    const satellite::State rk4_state =
        satellite::RungeKutta4Integrator{}.integrate({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
                                                    1.0, 0.5,
                                                    constant_velocity_derivative);
    if (!close(rk4_state.position.x, 1.0, tolerance) ||
        !close(rk4_state.position.y, 0.0, tolerance) ||
        !close(rk4_state.position.z, 0.0, tolerance) ||
        !close(rk4_state.velocity.x, 1.0, tolerance) ||
        !close(rk4_state.velocity.y, 0.0, tolerance) ||
        !close(rk4_state.velocity.z, 0.0, tolerance)) {
        std::cerr << "Runge-Kutta integrator did not respect the generic state derivative\n";
        return 1;
    }

    const double gravitational_parameter =
        satellite::constants::earth_gravitational_parameter_km3_s2;
    const double orbit_radius_km = 7000.0;
    const double circular_speed_km_s =
        std::sqrt(gravitational_parameter / orbit_radius_km);
    const double mean_motion_rad_s =
        std::sqrt(gravitational_parameter /
                  (orbit_radius_km * orbit_radius_km * orbit_radius_km));
    const double orbital_period_seconds =
        2.0 * satellite::constants::pi / mean_motion_rad_s;
    const double comparison_step_seconds = 60.0;
    const satellite::State circular_initial_state{
        {orbit_radius_km, 0.0, 0.0}, {0.0, circular_speed_km_s, 0.0}};
    const satellite::State analytical_final_state{
        {orbit_radius_km * std::cos(mean_motion_rad_s * orbital_period_seconds),
         orbit_radius_km * std::sin(mean_motion_rad_s * orbital_period_seconds),
         0.0},
        {-circular_speed_km_s *
             std::sin(mean_motion_rad_s * orbital_period_seconds),
         circular_speed_km_s *
             std::cos(mean_motion_rad_s * orbital_period_seconds),
         0.0}};
    const satellite::State euler_orbit_state = propagator.propagate_euler(
        circular_initial_state, orbital_period_seconds, comparison_step_seconds);
    const satellite::State rk4_orbit_state = propagator.propagate_rk4(
        circular_initial_state, orbital_period_seconds, comparison_step_seconds);

    constexpr double known_orbital_period_seconds = 5828.516638;
    constexpr double validation_interval_seconds = 60.0;
    constexpr double validation_step_seconds = 10.0;
    satellite::State sampled_orbit_state = circular_initial_state;
    double elapsed_seconds = 0.0;
    double maximum_radius_error_km = 0.0;
    double maximum_speed_error_km_s = 0.0;
    while (elapsed_seconds < orbital_period_seconds) {
        const double interval = std::min(validation_interval_seconds,
                                         orbital_period_seconds - elapsed_seconds);
        sampled_orbit_state = propagator.propagate_rk4(
            sampled_orbit_state, interval, validation_step_seconds);
        elapsed_seconds += interval;
        maximum_radius_error_km = std::max(
            maximum_radius_error_km,
            std::abs(magnitude(sampled_orbit_state.position) - orbit_radius_km));
        maximum_speed_error_km_s = std::max(
            maximum_speed_error_km_s,
            std::abs(magnitude(sampled_orbit_state.velocity) - circular_speed_km_s));
    }
    const double return_position_error_km =
        magnitude(sampled_orbit_state.position - circular_initial_state.position);

    if (!close(orbital_period_seconds, known_orbital_period_seconds, 0.001) ||
        maximum_radius_error_km > 0.001 ||
        maximum_speed_error_km_s > 0.000001 ||
        return_position_error_km > 0.001) {
        std::cerr << "Circular orbit validation failed: period="
                  << orbital_period_seconds << " s, max radius error="
                  << maximum_radius_error_km << " km, max speed error="
                  << maximum_speed_error_km_s << " km/s, return position error="
                  << return_position_error_km << " km\n";
        return 1;
    }

    const auto measure_errors = [&](const satellite::State& state) {
        const double position_error_km =
            magnitude(state.position - analytical_final_state.position);
        const double velocity_error_km_s =
            magnitude(state.velocity - analytical_final_state.velocity);
        const double energy_error_percent =
            100.0 * std::abs(specific_energy(state, gravitational_parameter) -
                             specific_energy(circular_initial_state,
                                             gravitational_parameter)) /
            std::abs(specific_energy(circular_initial_state,
                                     gravitational_parameter));
        const double angular_momentum_error_percent =
            100.0 * std::abs(angular_momentum_magnitude(state) -
                             angular_momentum_magnitude(circular_initial_state)) /
            angular_momentum_magnitude(circular_initial_state);
        return AccuracyErrors{position_error_km, velocity_error_km_s,
                      energy_error_percent,
                      angular_momentum_error_percent};
    };
    const auto euler_errors = measure_errors(euler_orbit_state);
    const auto rk4_errors = measure_errors(rk4_orbit_state);

    if (!(rk4_errors.position_error_km < euler_errors.position_error_km &&
          rk4_errors.velocity_error_km_s < euler_errors.velocity_error_km_s &&
          rk4_errors.energy_error_percent < euler_errors.energy_error_percent &&
          rk4_errors.angular_momentum_error_percent <
              euler_errors.angular_momentum_error_percent)) {
        std::cerr << "RK4 did not improve all one-orbit accuracy metrics over Euler\n";
        return 1;
    }

    std::cout << std::fixed << std::setprecision(6)
              << "One-orbit accuracy comparison (step: "
              << comparison_step_seconds << " s, period: "
              << orbital_period_seconds << " s)\n"
              << "Integrator | Position error (km) | Velocity error (km/s) | "
                 "Energy error (%) | Angular momentum error (%)\n"
              << "Euler      | " << euler_errors.position_error_km << " | "
              << euler_errors.velocity_error_km_s << " | "
              << euler_errors.energy_error_percent << " | "
              << euler_errors.angular_momentum_error_percent << "\n"
              << "RK4        | " << rk4_errors.position_error_km << " | "
              << rk4_errors.velocity_error_km_s << " | "
              << rk4_errors.energy_error_percent << " | "
              << rk4_errors.angular_momentum_error_percent << "\n";

    return 0;
}
