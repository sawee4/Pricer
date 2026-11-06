#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <tuple>

#include "pricer/analytic_pricer.hpp"
#include "pricer/black_scholes.hpp"

using namespace pricer;
using namespace pricer::black_scholes;

namespace {

double price(OptionType type, const Inputs& in) { return vanillaPrice(type, in); }

// Central finite differences. The spot step is scaled by S*sigma*sqrt(T), the natural
// width of the price curve: a fixed step would be too coarse for short-dated, low-vol
// options (large truncation error O(h^2)) and too fine for long-dated ones (round-off
// error O(eps/h) for delta, O(eps/h^2) for gamma).
constexpr double kDeltaStep = 1e-4;   // in units of S sigma sqrt(T)
constexpr double kGammaStep = 1e-3;   // larger: round-off of a 2nd difference grows as 1/h^2
constexpr double kParamStep = 1e-5;   // absolute step for sigma, T and r

double spotScale(const Inputs& in) {
    return in.spot * in.volatility * std::sqrt(in.maturity);
}

double fdDelta(OptionType type, Inputs in) {
    const double h = kDeltaStep * spotScale(in);
    Inputs up = in, down = in;
    up.spot += h;
    down.spot -= h;
    return (price(type, up) - price(type, down)) / (2.0 * h);
}

double fdGamma(OptionType type, Inputs in) {
    const double h = kGammaStep * spotScale(in);
    Inputs up = in, down = in;
    up.spot += h;
    down.spot -= h;
    return (price(type, up) - 2.0 * price(type, in) + price(type, down)) / (h * h);
}

double fdVega(OptionType type, Inputs in) {
    Inputs up = in, down = in;
    up.volatility += kParamStep;
    down.volatility -= kParamStep;
    return (price(type, up) - price(type, down)) / (2.0 * kParamStep);
}

double fdTheta(OptionType type, Inputs in) {
    // theta = dV/dt = -dV/dT (time passing shortens the maturity).
    Inputs up = in, down = in;
    up.maturity += kParamStep;
    down.maturity -= kParamStep;
    return -(price(type, up) - price(type, down)) / (2.0 * kParamStep);
}

double fdRho(OptionType type, Inputs in) {
    Inputs up = in, down = in;
    up.rate += kParamStep;
    down.rate -= kParamStep;
    return (price(type, up) - price(type, down)) / (2.0 * kParamStep);
}

double tolerance(double value, double rel) { return rel * std::max(1.0, std::abs(value)); }

}  // namespace

// ---------------------------------------------------------------------------
// Textbook values: Hull, Examples 19.1-19.6 (S=49, K=50, r=5%, sigma=20%, T=20 weeks)
// ---------------------------------------------------------------------------
TEST(Greeks, HullTextbookValues) {
    const Inputs in{49.0, 50.0, 0.05, 0.0, 0.20, 20.0 / 52.0};
    const Greeks g = vanillaGreeks(OptionType::Call, in);
    EXPECT_NEAR(g.delta, 0.522, 5e-4);
    EXPECT_NEAR(g.gamma, 0.066, 5e-4);
    EXPECT_NEAR(g.vega, 12.1, 0.05);
    EXPECT_NEAR(g.theta, -4.31, 0.005);
    EXPECT_NEAR(g.rho, 8.91, 0.005);
}

// ---------------------------------------------------------------------------
// Analytic Greeks vs central finite differences, over a parameter grid
// ---------------------------------------------------------------------------
class GreeksVsFiniteDifferences
    : public ::testing::TestWithParam<std::tuple<OptionType, double, double, double, double>> {};

TEST_P(GreeksVsFiniteDifferences, Agree) {
    const auto [type, spot, dividend, vol, maturity] = GetParam();
    const Inputs in{spot, 100.0, 0.05, dividend, vol, maturity};
    const Greeks g = vanillaGreeks(type, in);

    EXPECT_NEAR(g.delta, fdDelta(type, in), 1e-8);
    EXPECT_NEAR(g.gamma, fdGamma(type, in), tolerance(g.gamma, 1e-6));
    EXPECT_NEAR(g.vega, fdVega(type, in), tolerance(g.vega, 1e-6));
    EXPECT_NEAR(g.theta, fdTheta(type, in), tolerance(g.theta, 1e-6));
    EXPECT_NEAR(g.rho, fdRho(type, in), tolerance(g.rho, 1e-6));
}

INSTANTIATE_TEST_SUITE_P(
    Grid, GreeksVsFiniteDifferences,
    ::testing::Combine(::testing::Values(OptionType::Call, OptionType::Put),
                       ::testing::Values(70.0, 95.0, 100.0, 105.0, 140.0),  // spot (K = 100)
                       ::testing::Values(0.0, 0.03),                        // dividend
                       ::testing::Values(0.1, 0.3, 0.6),                    // vol
                       ::testing::Values(0.1, 1.0, 5.0)));                  // maturity

// ---------------------------------------------------------------------------
// Model-consistency checks
// ---------------------------------------------------------------------------
TEST(Greeks, SatisfyBlackScholesPde) {
    // theta + (r - q) S delta + 1/2 sigma^2 S^2 gamma - r V = 0
    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        for (double spot : {60.0, 100.0, 150.0}) {
            const Inputs in{spot, 100.0, 0.04, 0.02, 0.25, 1.5};
            const Greeks g = vanillaGreeks(type, in);
            const double v = price(type, in);
            const double pde = g.theta + (in.rate - in.dividend) * spot * g.delta +
                               0.5 * in.volatility * in.volatility * spot * spot * g.gamma -
                               in.rate * v;
            EXPECT_NEAR(pde, 0.0, 1e-12);
        }
    }
}

TEST(Greeks, CallPutRelationsFromParity) {
    // Differentiating C - P = S e^{-qT} - K e^{-rT}.
    const Inputs in{105.0, 100.0, 0.03, 0.01, 0.3, 2.0};
    const Greeks c = vanillaGreeks(OptionType::Call, in);
    const Greeks p = vanillaGreeks(OptionType::Put, in);
    EXPECT_NEAR(c.delta - p.delta, std::exp(-0.01 * 2.0), 1e-14);
    EXPECT_NEAR(c.gamma, p.gamma, 1e-15);
    EXPECT_NEAR(c.vega, p.vega, 1e-12);
    EXPECT_NEAR(c.rho - p.rho, 100.0 * 2.0 * std::exp(-0.03 * 2.0), 1e-12);
    EXPECT_NEAR(c.theta - p.theta,
                0.01 * 105.0 * std::exp(-0.01 * 2.0) - 0.03 * 100.0 * std::exp(-0.03 * 2.0),
                1e-12);
}

TEST(Greeks, SignsAndBounds) {
    const Inputs in{100.0, 100.0, 0.05, 0.0, 0.2, 1.0};
    const Greeks c = vanillaGreeks(OptionType::Call, in);
    const Greeks p = vanillaGreeks(OptionType::Put, in);
    EXPECT_GT(c.delta, 0.0);
    EXPECT_LT(c.delta, 1.0);
    EXPECT_LT(p.delta, 0.0);
    EXPECT_GT(p.delta, -1.0);
    EXPECT_GT(c.gamma, 0.0);
    EXPECT_GT(c.vega, 0.0);
    EXPECT_LT(c.theta, 0.0);  // no dividend: a call loses value as time passes
    EXPECT_GT(c.rho, 0.0);
    EXPECT_LT(p.rho, 0.0);
}

TEST(Greeks, DeterministicLimit) {
    // sigma = 0, deep in the money call: V = S e^{-qT} - K e^{-rT}.
    const Inputs in{120.0, 100.0, 0.05, 0.02, 0.0, 1.0};
    const Greeks g = vanillaGreeks(OptionType::Call, in);
    EXPECT_DOUBLE_EQ(g.delta, std::exp(-0.02));
    EXPECT_DOUBLE_EQ(g.gamma, 0.0);
    EXPECT_DOUBLE_EQ(g.vega, 0.0);
    EXPECT_NEAR(g.theta, fdTheta(OptionType::Call, in), 1e-6);
    EXPECT_NEAR(g.rho, fdRho(OptionType::Call, in), 1e-6);

    const Greeks otmPut = vanillaGreeks(OptionType::Put, in);
    EXPECT_DOUBLE_EQ(otmPut.delta, 0.0);
    EXPECT_DOUBLE_EQ(otmPut.theta, 0.0);
    EXPECT_DOUBLE_EQ(otmPut.rho, 0.0);
}

TEST(Greeks, AtTheMoneyGammaExplodesAsMaturityShrinks) {
    // gamma_ATM ~ n(0) / (S sigma sqrt T): the pin risk of short-dated options.
    double previous = 0.0;
    for (double t : {1.0, 0.1, 0.01, 0.001}) {
        const double gamma = vanillaGreeks(OptionType::Call, {100.0, 100.0, 0.0, 0.0, 0.2, t}).gamma;
        EXPECT_GT(gamma, previous);
        previous = gamma;
    }
}

// ---------------------------------------------------------------------------
// AnalyticPricer::greeks
// ---------------------------------------------------------------------------
TEST(Greeks, AnalyticPricerDispatch) {
    const AnalyticPricer pricer;
    const MarketData market(49.0, 0.05, 0.0, 0.20);
    const Option call(std::make_unique<CallPayoff>(50.0), 20.0 / 52.0);
    EXPECT_NEAR(pricer.greeks(call, market).delta, 0.522, 5e-4);

    const Option digital(std::make_unique<DigitalCallPayoff>(50.0), 1.0);
    EXPECT_THROW((void)pricer.greeks(digital, market), std::invalid_argument);
    const Option american(std::make_unique<PutPayoff>(50.0), 1.0, ExerciseStyle::American);
    EXPECT_THROW((void)pricer.greeks(american, market), std::invalid_argument);
}
