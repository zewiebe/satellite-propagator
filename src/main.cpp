#include "satellite/propagator.h"
#include "satellite/constants.h"

#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>

namespace {
constexpr double default_duration_seconds = 3600.0;
constexpr double default_step_seconds = 10.0;

void print_usage(const char* program) {
    std::cerr << "Usage: " << program << " [duration_seconds] [step_seconds]\n";
}
}

int main(int argc, char* argv[]) {
    try {
        const double duration = argc > 1 ? std::stod(argv[1]) : default_duration_seconds;
        const double step = argc > 2 ? std::stod(argv[2]) : default_step_seconds;
        if (argc > 3) {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }

        const satellite::State initial_state{
            {7000.0, 0.0, 0.0},
            {0.0, 7.54605329, 0.0}
        };
        satellite::Propagator propagator(
            satellite::constants::earth_gravitational_parameter_km3_s2);
        const satellite::State final_state = propagator.propagate(initial_state, duration, step);

        std::cout << std::fixed << std::setprecision(6)
                  << "duration_seconds," << duration << "\n"
                  << "step_seconds," << step << "\n"
                  << "final_position_km," << final_state.position.x << ','
                  << final_state.position.y << ',' << final_state.position.z << "\n"
                  << "final_velocity_km_s," << final_state.velocity.x << ','
                  << final_state.velocity.y << ',' << final_state.velocity.z << "\n";
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
