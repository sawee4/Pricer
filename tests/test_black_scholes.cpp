#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <tuple>

#include "pricer/black_scholes.hpp"
#include "pricer/normal.hpp"

using namespace pricer;
using namespace pricer::black_scholes;

namespace {

// Hull, "Options, Futures and Other Derivatives", Example 15.6.
constexpr Inputs kHull{42.0, 40.0, 0.10, 0.0, 0.20, 0.5};

}  // namespace

TEST(BlackScholes, HullTextbookExample) {
    // Textbook values (rounded to the cent): call 4.76, put 0.81.
    EXPECT_NEAR(vanillaPrice(OptionType::Call, kHull), 4.76, 0.005);
    EXPECT_NEAR(vanillaPrice(OptionType::Put, kHull), 0.81, 0.005);
    // Full precision regression values.
    EXPECT_NEAR(vanillaPrice(OptionType::Call, kHull), 4.759422392871536, 1e-12);
    EXPECT_NEAR(vanillaPrice(OptionType::Put, kHull), 0.808599372900093, 1e-12);
}

TEST(BlackScholes, HullD1D2) {
    const auto [d1, d2] = computeD1D2(kHull);
    EXPECT_NEAR(d1, 0.7693, 1e-4);
    EXPECT_NEAR(d2, 0.6278, 1e-4);
}

TEST(BlackScholes, HullIndexOptionWithDividendYield) {
    // Hull, index option example: S=930, K=900, r=8%, q=3%, sigma=20%, T=2 months -> 51.83.
    const Inputs in{930.0, 900.0, 0.08, 0.03, 0.20, 2.0 / 12.0};
    EXPECT_NEAR(vanillaPrice(OptionType::Call, in), 51.83, 0.005);
}

// ---------------------------------------------------------------------------
// Put-call parity: C - P = S e^{-qT} - K e^{-rT}, over a grid of parameters.
// ---------------------------------------------------------------------------
class PutCallParity
    : public ::testing::TestWithParam<std::tuple<double, double, double, double, double>> {};

TEST_P(PutCallParity, Holds) {
    const auto [spot, strike, rate, dividend, vol] = GetParam();
    for (double maturity : {0.01, 0.5, 2.0, 10.0}) {
        const Inputs in{spot, strike, rate, dividend, vol, maturity};
        const double lhs = vanillaPrice(OptionType::Call, in) - vanillaPrice(OptionType::Put, in);
        const double rhs = spot * std::exp(-dividend * maturity) - strike * std::exp(-rate * maturity);
        EXPECT_NEAR(lhs, rhs, 1e-10 * std::max(1.0, spot));
    }
}

INSTANTIATE_TEST_SUITE_P(
    Grid, PutCallParity,
    ::testing::Combine(::testing::Values(50.0, 100.0, 150.0),   // spot
                       ::testing::Values(80.0, 100.0, 120.0),   // strike
                       ::testing::Values(-0.01, 0.0, 0.05),     // rate
                       ::testing::Values(0.0, 0.03),            // dividend
                       ::testing::Values(0.05, 0.2, 0.8)));     // vol

// ---------------------------------------------------------------------------
// No-arbitrage bounds and monotonicity
// ---------------------------------------------------------------------------
TEST(BlackScholes, CallWithinNoArbitrageBounds) {
    for (double spot : {60.0, 100.0, 140.0}) {
        const Inputs in{spot, 100.0, 0.05, 0.02, 0.3, 1.0};
        const double call = vanillaPrice(OptionType::Call, in);
        const double lower = std::max(spot * std::exp(-0.02) - 100.0 * std::exp(-0.05), 0.0);
        EXPECT_GE(call, lower);
        EXPECT_LE(call, spot * std::exp(-0.02));
    }
}

TEST(BlackScholes, PriceIncreasesWithVolatility) {
    double previous = 0.0;
    for (double vol : {0.05, 0.1, 0.2, 0.4, 0.8}) {
        const double call = vanillaPrice(OptionType::Call, {100.0, 100.0, 0.05, 0.0, vol, 1.0});
        EXPECT_GT(call, previous);
        previous = call;
    }
}

// ---------------------------------------------------------------------------
// Limit cases: sigma -> 0 and T -> 0
// ---------------------------------------------------------------------------
TEST(BlackScholes, ZeroVolatilityIsDiscountedForwardIntrinsic) {
    const Inputs in{100.0, 100.0, 0.05, 0.01, 0.0, 1.0};
    const double df = std::exp(-0.05);
    const double fwd = 100.0 * std::exp(0.04);
    EXPECT_DOUBLE_EQ(vanillaPrice(OptionType::Call, in), df * (fwd - 100.0));
    EXPECT_DOUBLE_EQ(vanillaPrice(OptionType::Put, in), 0.0);
}

TEST(BlackScholes, VanishingVolatilityConvergesContinuously) {
    // Formula path (sigma small but above threshold) must agree with the sigma=0 branch.
    for (double strike : {90.0, 100.0, 110.0}) {
        const Inputs zero{100.0, strike, 0.05, 0.01, 0.0, 1.0};
        const Inputs tiny{100.0, strike, 0.05, 0.01, 1e-8, 1.0};
        EXPECT_NEAR(vanillaPrice(OptionType::Call, tiny), vanillaPrice(OptionType::Call, zero),
                    1e-7);
        EXPECT_NEAR(vanillaPrice(OptionType::Put, tiny), vanillaPrice(OptionType::Put, zero),
                    1e-7);
    }
}

TEST(BlackScholes, ZeroMaturityIsIntrinsicValue) {
    EXPECT_DOUBLE_EQ(vanillaPrice(OptionType::Call, {110.0, 100.0, 0.05, 0.0, 0.2, 0.0}), 10.0);
    EXPECT_DOUBLE_EQ(vanillaPrice(OptionType::Put, {110.0, 100.0, 0.05, 0.0, 0.2, 0.0}), 0.0);
    EXPECT_DOUBLE_EQ(vanillaPrice(OptionType::Put, {90.0, 100.0, 0.05, 0.0, 0.2, 0.0}), 10.0);
}

TEST(BlackScholes, ShortMaturityConvergesToIntrinsic) {
    // ATM the time value decays like sigma*S*sqrt(T/(2 pi)); away from the money it is
    // exponentially small.
    const double t = 1e-8;
    EXPECT_NEAR(vanillaPrice(OptionType::Call, {110.0, 100.0, 0.05, 0.0, 0.2, t}), 10.0, 1e-6);
    const double atm = vanillaPrice(OptionType::Call, {100.0, 100.0, 0.0, 0.0, 0.2, t});
    EXPECT_NEAR(atm, 0.2 * 100.0 * std::sqrt(t) * pricer::math::kInvSqrt2Pi, 1e-9);
}

TEST(BlackScholes, DeepOutOfTheMoneyIsTinyButFiniteAndNonNegative) {
    const Inputs farOtm{1.0, 1000.0, 0.05, 0.0, 0.2, 1.0};
    const double call = vanillaPrice(OptionType::Call, farOtm);
    EXPECT_TRUE(std::isfinite(call));
    EXPECT_GE(call, 0.0);
    EXPECT_LT(call, 1e-100);
}

// ---------------------------------------------------------------------------
// Digitals
// ---------------------------------------------------------------------------
TEST(BlackScholes, DigitalCallPlusDigitalPutIsZeroCouponBond) {
    const Inputs in{100.0, 95.0, 0.05, 0.02, 0.25, 1.5};
    EXPECT_NEAR(digitalPrice(OptionType::Call, in) + digitalPrice(OptionType::Put, in),
                std::exp(-0.05 * 1.5), 1e-14);
}

TEST(BlackScholes, DigitalCallIsMinusStrikeDerivativeOfCall) {
    // Breeden-Litzenberger: digital call = -dC/dK. Central difference, O(h^2) error.
    const double h = 1e-4;
    Inputs up{100.0, 100.0 + h, 0.05, 0.02, 0.25, 1.0};
    Inputs down{100.0, 100.0 - h, 0.05, 0.02, 0.25, 1.0};
    const double fd =
        -(vanillaPrice(OptionType::Call, up) - vanillaPrice(OptionType::Call, down)) / (2.0 * h);
    EXPECT_NEAR(digitalPrice(OptionType::Call, {100.0, 100.0, 0.05, 0.02, 0.25, 1.0}), fd, 1e-7);
}

TEST(BlackScholes, RejectsInvalidInputs) {
    EXPECT_THROW((void)vanillaPrice(OptionType::Call, {-1.0, 100.0, 0.05, 0.0, 0.2, 1.0}),
                 std::invalid_argument);
    EXPECT_THROW((void)vanillaPrice(OptionType::Call, {100.0, 100.0, 0.05, 0.0, -0.2, 1.0}),
                 std::invalid_argument);
    EXPECT_THROW((void)computeD1D2({100.0, 100.0, 0.05, 0.0, 0.0, 1.0}), std::domain_error);
}
