#pragma once

namespace satellite::constants {

// Distance units are kilometers, time units are seconds, and mass units are kilograms.
inline constexpr double earth_radius_km = 6378.137;
inline constexpr double earth_gravitational_parameter_km3_s2 = 398600.4418;

inline constexpr double pi = 3.141592653589793238462643383279502884;
inline constexpr double degrees_to_radians = pi / 180.0;
inline constexpr double radians_to_degrees = 180.0 / pi;

} // namespace satellite::constants
