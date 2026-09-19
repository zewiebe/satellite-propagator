#pragma once

namespace satellite {

struct Vector3 {
    double x;
    double y;
    double z;

    Vector3 operator+(const Vector3& other) const;
    Vector3 operator-(const Vector3& other) const;
    Vector3 operator*(double scalar) const;
};

struct State {
    Vector3 position;
    Vector3 velocity;
};

class Propagator {
public:
    explicit Propagator(double gravitational_parameter_km3_s2);

    State propagate(const State& initial_state, double duration_seconds,
                    double step_seconds) const;

    Vector3 acceleration(const Vector3& position) const;

private:
    double gravitational_parameter_;
};

}
