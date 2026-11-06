include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2
    GIT_SHALLOW    TRUE
)
# On MSVC, link GoogleTest against the same (dynamic) CRT as our code.
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(pricer_tests
    test_payoff.cpp
    test_market_data.cpp
    test_option.cpp
    test_normal.cpp
    test_black_scholes.cpp
    test_analytic_pricer.cpp
    test_greeks.cpp
    test_random.cpp
    test_monte_carlo.cpp
)
target_link_libraries(pricer_tests PRIVATE pricer::pricer pricer_warnings GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(pricer_tests DISCOVERY_MODE PRE_TEST)
