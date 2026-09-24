#pragma once

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

class Propagator {
public:
    // gravitational_parameter_km3_s2 is mu in km^3/s^2.
    explicit Propagator(double gravitational_parameter_km3_s2);

    // duration_seconds and step_seconds are in seconds (s).
    State propagate(const State& initial_state, double duration_seconds,
                    double step_seconds) const;

    // Returns Earth-centered two-body acceleration in km/s^2. The supplied
    // gravitational parameter is the only force parameter; no perturbations
    // are modeled.
    Vector3 acceleration(const Vector3& position) const;

private:
    double gravitational_parameter_;
};

}
