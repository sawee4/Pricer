cmake_minimum_required(VERSION 3.16)

project(pricer
    VERSION 0.1.0
    DESCRIPTION "Option pricing library: analytic, Monte Carlo and binomial pricers"
    LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(PRICER_BUILD_TESTS          "Build the unit tests (fetches GoogleTest)" ON)
option(PRICER_BUILD_EXAMPLES       "Build the example programs"                ON)
option(PRICER_WARNINGS_AS_ERRORS   "Treat compiler warnings as errors"         OFF)

# Default to Release for single-config generators (Makefiles, Ninja).
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# Warnings, applied to our own targets only (not to GoogleTest)
# ---------------------------------------------------------------------------
add_library(pricer_warnings INTERFACE)
if(MSVC)
    target_compile_options(pricer_warnings INTERFACE /W4 /permissive- /utf-8)
    if(PRICER_WARNINGS_AS_ERRORS)
        target_compile_options(pricer_warnings INTERFACE /WX)
    endif()
else()
    target_compile_options(pricer_warnings INTERFACE
        -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
        -Wold-style-cast -Wnon-virtual-dtor -Woverloaded-virtual)
    if(PRICER_WARNINGS_AS_ERRORS)
        target_compile_options(pricer_warnings INTERFACE -Werror)
    endif()
endif()

# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------
add_library(pricer
    src/payoff.cpp
    src/black_scholes.cpp
    src/analytic_pricer.cpp
    src/normal.cpp
    src/random.cpp
    src/monte_carlo_pricer.cpp
    src/binomial_pricer.cpp
    src/monte_carlo_greeks.cpp
    src/implied_vol.cpp
    src/path_payoff.cpp
    src/path_dependent_pricer.cpp
    src/exotic_formulas.cpp
    src/market_data.cpp
    src/option.cpp
)
add_library(pricer::pricer ALIAS pricer)

target_include_directories(pricer PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_compile_features(pricer PUBLIC cxx_std_17)
target_link_libraries(pricer PRIVATE pricer_warnings)

# No floating-point contraction: GCC fuses a * b + c into one FMA instruction whenever the
# target has FMA (it does so even in ISO C++ mode), which skips the rounding of a * b. Results
# would then depend on the CPU flags the compiler targets (bit-reproducibility across machines
# is lost, and Welford's variance of identical samples is no longer exactly zero). PUBLIC
# because the statistics accumulators are inline in the headers. MSVC does not contract by
# default (/fp:precise).
if(NOT MSVC)
    target_compile_options(pricer PUBLIC -ffp-contract=off)
endif()

# ---------------------------------------------------------------------------
# Tests and examples
# ---------------------------------------------------------------------------
if(PRICER_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(PRICER_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
