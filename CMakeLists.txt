cmake_minimum_required(VERSION 3.16)

project(satellite-propagator VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(propagator
    src/propagator.cpp
)

target_include_directories(propagator PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(propagator PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive->
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
)

add_executable(satellite-propagator
    src/main.cpp
)

target_link_libraries(satellite-propagator PRIVATE propagator)

enable_testing()

include(CTest)

if(BUILD_TESTING)
    add_executable(propagator-tests
        tests/propagator_tests.cpp
    )
    target_link_libraries(propagator-tests PRIVATE propagator)
    add_test(NAME propagator-tests COMMAND propagator-tests)
endif()
