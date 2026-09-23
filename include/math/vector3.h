#pragma once

#include <cmath>

namespace math {

struct Vector3 {
    double x;
    double y;
    double z;

    constexpr Vector3() noexcept : x(0.0), y(0.0), z(0.0) {}
    constexpr Vector3(double x_value, double y_value, double z_value) noexcept
        : x(x_value), y(y_value), z(z_value) {}

    // Vector addition: c = a + b
    // Each component is summed independently:
    // c.x = a.x + b.x, c.y = a.y + b.y, c.z = a.z + b.z
    constexpr Vector3 operator+(const Vector3& rhs) const noexcept {
        return {x + rhs.x, y + rhs.y, z + rhs.z};
    }

    // Vector subtraction: c = a - b
    // Each component is differenced independently:
    // c.x = a.x - b.x, c.y = a.y - b.y, c.z = a.z - b.z
    constexpr Vector3 operator-(const Vector3& rhs) const noexcept {
        return {x - rhs.x, y - rhs.y, z - rhs.z};
    }

    // Vector negation: -v = (-v.x, -v.y, -v.z)
    // This flips the direction of the vector while preserving magnitude.
    constexpr Vector3 operator-() const noexcept {
        return {-x, -y, -z};
    }

    // In-place vector addition: this = this + rhs
    constexpr Vector3& operator+=(const Vector3& rhs) noexcept {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    // In-place vector subtraction: this = this - rhs
    constexpr Vector3& operator-=(const Vector3& rhs) noexcept {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    // Scalar multiplication: v * s = (v.x*s, v.y*s, v.z*s)
    // Scales the vector's length by s while keeping direction.
    constexpr Vector3 operator*(double scalar) const noexcept {
        return {x * scalar, y * scalar, z * scalar};
    }

    // Scalar division: v / s = (v.x/s, v.y/s, v.z/s)
    // Divides each component by a scalar value.
    constexpr Vector3 operator/(double scalar) const noexcept {
        return {x / scalar, y / scalar, z / scalar};
    }

    // In-place scalar multiplication: this = this * scalar
    constexpr Vector3& operator*=(double scalar) noexcept {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    // In-place scalar division: this = this / scalar
    constexpr Vector3& operator/=(double scalar) noexcept {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }

    // Euclidean norm: ||v|| = sqrt(v.x^2 + v.y^2 + v.z^2)
    // Measures the vector length from the origin.
    double magnitude() const noexcept {
        return std::sqrt(x * x + y * y + z * z);
    }

    // Squared magnitude: ||v||^2 = x^2 + y^2 + z^2
    // Used to avoid a sqrt when only relative magnitude is needed.
    constexpr double magnitude_squared() const noexcept {
        return x * x + y * y + z * z;
    }

    // Unit vector in the same direction as this vector:
    // v_hat = v / ||v|| when ||v|| != 0; otherwise returns zero vector.
    Vector3 normalized() const noexcept {
        const double length = magnitude();
        return length == 0.0 ? Vector3{} : *this / length;
    }

    // Converts this vector to unit length in-place.
    void normalize() noexcept {
        const double length = magnitude();
        if (length != 0.0) {
            *this /= length;
        }
    }

    // Dot product: a · b = a.x*b.x + a.y*b.y + a.z*b.z
    // Measures alignment between two vectors.
    constexpr double dot(const Vector3& rhs) const noexcept {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    // Cross product: a × b
    // Produces a vector perpendicular to both a and b:
    // (a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x)
    constexpr Vector3 cross(const Vector3& rhs) const noexcept {
        return {y * rhs.z - z * rhs.y,
                z * rhs.x - x * rhs.z,
                x * rhs.y - y * rhs.x};
    }
};

 // Scalar-vector multiplication for commutative form:
 // s * v = v * s
 constexpr Vector3 operator*(double scalar, const Vector3& vector) noexcept {
     return vector * scalar;
 }

} // namespace math