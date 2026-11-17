#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "pricer/implied_vol.hpp"

using namespace pricer;
using black_scholes::Inputs;

TEST(ImpliedVol, HullExample) {
    // Hull: European call, S=21, K=20, r=10%, T=3 months, price 1.875 -> sigma ~ 23.5%.
    const ImpliedVolResult r =
        impliedVolatility(OptionType::Call, 1.875, {21.0, 20.0, 0.10, 0.0, 0.0, 0.25});
    EXPECT_NEAR(r.volatility, 0.235, 5e-4);
    EXPECT_LE(r.iterations, 8U);  // Newton: quadratic convergence
}

TEST(ImpliedVol, RoundTripOverWideGrid) {
    // Price with a known sigma, invert, compare. Skipped: cases whose time value is below
    // ~1e-10 relative, where the price is numerically equal to its sigma -> 0 bound and
    // carries no information about sigma (any solver must reject them).
    int checked = 0;
    for (double sigma : {0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 3.0}) {
        for (double strike : {50.0, 80.0, 100.0, 120.0, 200.0}) {
            for (double maturity : {0.02, 0.5, 5.0}) {
                for (OptionType type : {OptionType::Call, OptionType::Put}) {
                    const Inputs in{100.0, strike, 0.03, 0.01, sigma, maturity};
                    const double price = black_scholes::vanillaPrice(type, in);
                    Inputs zeroVol = in;
                    zeroVol.volatility = 0.0;
                    const double timeValue = price - black_scholes::vanillaPrice(type, zeroVol);
                    if (timeValue < 1e-10 * std::max(1.0, price)) continue;

                    const ImpliedVolResult r = impliedVolatility(type, price, in);
                    // Accuracy in sigma ~ price tolerance / vega.
                    const double vega = black_scholes::vanillaGreeks(type, in).vega;
                    const double tol = std::max(1e-10, 1e-11 * std::max(1.0, price) / vega);
                    EXPECT_NEAR(r.volatility, sigma, tol)
                        << "sigma=" << sigma << " K=" << strike << " T=" << maturity;
                    ++checked;
                }
            }
        }
    }
    EXPECT_GT(checked, 150);
}

TEST(ImpliedVol, CallAndPutGiveSameVolatility) {
    // Put-call parity holds for any sigma, so a parity-consistent call/put pair has one vol.
    const Inputs in{100.0, 110.0, 0.04, 0.02, 0.37, 1.5};
    const double call = black_scholes::vanillaPrice(OptionType::Call, in);
    const double put = call - 100.0 * std::exp(-0.02 * 1.5) + 110.0 * std::exp(-0.04 * 1.5);
    EXPECT_NEAR(impliedVolatility(OptionType::Call, call, in).volatility,
                impliedVolatility(OptionType::Put, put, in).volatility, 1e-10);
}

TEST(ImpliedVol, RejectsArbitrageablePrices) {
    const Inputs in{100.0, 100.0, 0.05, 0.0, 0.0, 1.0};
    const double lower = 100.0 - 100.0 * std::exp(-0.05);
    EXPECT_THROW((void)impliedVolatility(OptionType::Call, lower - 0.01, in), std::domain_error);
    EXPECT_THROW((void)impliedVolatility(OptionType::Call, 100.0, in), std::domain_error);
    EXPECT_THROW((void)impliedVolatility(OptionType::Put, -1.0, in), std::domain_error);
    const Inputs expired{100.0, 100.0, 0.05, 0.0, 0.0, 0.0};
    EXPECT_THROW((void)impliedVolatility(OptionType::Call, 5.0, expired), std::invalid_argument);
}

TEST(ImpliedVol, RespectsSearchBounds) {
    // True vol 3.0 but the caller caps the search at 2.0: reported, not silently clamped.
    const Inputs in{100.0, 100.0, 0.0, 0.0, 3.0, 1.0};
    const double price = black_scholes::vanillaPrice(OptionType::Call, in);
    ImpliedVolSettings capped;
    capped.upperBound = 2.0;
    EXPECT_THROW((void)impliedVolatility(OptionType::Call, price, in, capped), std::domain_error);
}
