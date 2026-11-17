#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "pricer/analytic_pricer.hpp"
#include "pricer/binomial_pricer.hpp"

using namespace pricer;

namespace {

double tree(std::size_t steps, const Option& option, const MarketData& market) {
    return BinomialPricer(steps).price(option, market).price;
}

double analytic(const Option& option, const MarketData& market) {
    return AnalyticPricer().price(option, market).price;
}

Option european(std::unique_ptr<Payoff> payoff, double maturity) {
    return Option(std::move(payoff), maturity, ExerciseStyle::European);
}

Option american(std::unique_ptr<Payoff> payoff, double maturity) {
    return Option(std::move(payoff), maturity, ExerciseStyle::American);
}

}  // namespace

// ---------------------------------------------------------------------------
// Textbook values: Hull, Example 21.1 and Table 21.1 (American put, S=K=50, r=10%,
// sigma=40%, T=5 months): 4.49 with 5 steps, then 4.263 / 4.272 / 4.278 / 4.283.
// ---------------------------------------------------------------------------
TEST(Binomial, HullAmericanPut) {
    const MarketData market(50.0, 0.10, 0.0, 0.40);
    const Option put = american(std::make_unique<PutPayoff>(50.0), 5.0 / 12.0);
    EXPECT_NEAR(tree(5, put, market), 4.49, 0.005);
    EXPECT_NEAR(tree(30, put, market), 4.263, 0.0005);
    EXPECT_NEAR(tree(50, put, market), 4.272, 0.0005);
    EXPECT_NEAR(tree(100, put, market), 4.278, 0.0005);
    EXPECT_NEAR(tree(500, put, market), 4.283, 0.0005);
}

// ---------------------------------------------------------------------------
// Convergence of the European tree to Black-Scholes
// ---------------------------------------------------------------------------
TEST(Binomial, EuropeanConvergesToBlackScholes) {
    const MarketData market(42.0, 0.10, 0.0, 0.20);
    for (bool isCall : {true, false}) {
        std::unique_ptr<Payoff> payoff;
        if (isCall) {
            payoff = std::make_unique<CallPayoff>(40.0);
        } else {
            payoff = std::make_unique<PutPayoff>(40.0);
        }
        const Option option = european(std::move(payoff), 0.5);
        const double bs = analytic(option, market);

        // The error oscillates with n, so compare the worst case over a window of n.
        double coarse = 0.0, fine = 0.0;
        for (std::size_t n = 50; n < 60; ++n) {
            coarse = std::max(coarse, std::abs(tree(n, option, market) - bs));
        }
        for (std::size_t n = 1000; n < 1010; ++n) {
            fine = std::max(fine, std::abs(tree(n, option, market) - bs));
        }
        EXPECT_LT(fine, 3e-3);
        EXPECT_LT(fine, coarse / 10.0);  // O(1/n): 20x more steps, at least 10x less error
    }
}

TEST(Binomial, AveragingConsecutiveStepsDampsOscillation) {
    // The CRR error has an odd/even component (averaging n and n+1 cancels most of it) and a
    // slower one driven by where ln(K) falls between terminal nodes (averaging does not fix
    // that one). Over a window of n, the worst averaged error is at least halved.
    const MarketData market(100.0, 0.05, 0.02, 0.25);
    const Option call = european(std::make_unique<CallPayoff>(105.0), 1.0);
    const double bs = analytic(call, market);
    double worstRaw = 0.0, worstAveraged = 0.0;
    double current = tree(200, call, market);
    for (std::size_t n = 200; n < 300; ++n) {
        const double next = tree(n + 1, call, market);
        worstRaw = std::max(worstRaw, std::abs(current - bs));
        worstAveraged = std::max(worstAveraged, std::abs(0.5 * (current + next) - bs));
        current = next;
    }
    EXPECT_LT(worstAveraged, 0.5 * worstRaw);
}

class EuropeanTreeMatchesBlackScholes
    : public ::testing::TestWithParam<std::tuple<double, double, double, double>> {};

TEST_P(EuropeanTreeMatchesBlackScholes, CloseToBlackScholesWith2000Steps) {
    const auto [strike, dividend, vol, maturity] = GetParam();
    const MarketData market(100.0, 0.04, dividend, vol);
    const Option call = european(std::make_unique<CallPayoff>(strike), maturity);
    const Option put = european(std::make_unique<PutPayoff>(strike), maturity);
    EXPECT_NEAR(tree(2000, call, market), analytic(call, market), 1e-2 * vol * std::sqrt(maturity));
    EXPECT_NEAR(tree(2000, put, market), analytic(put, market), 1e-2 * vol * std::sqrt(maturity));
}

INSTANTIATE_TEST_SUITE_P(Grid, EuropeanTreeMatchesBlackScholes,
                         ::testing::Combine(::testing::Values(80.0, 100.0, 120.0),
                                            ::testing::Values(0.0, 0.03),
                                            ::testing::Values(0.15, 0.4),
                                            ::testing::Values(0.25, 2.0)));

// ---------------------------------------------------------------------------
// Early exercise properties
// ---------------------------------------------------------------------------
TEST(Binomial, AmericanPutWorthAtLeastEuropeanPut) {
    for (double spot : {70.0, 90.0, 100.0, 110.0, 130.0}) {
        const MarketData market(spot, 0.06, 0.0, 0.3);
        const double us = tree(500, american(std::make_unique<PutPayoff>(100.0), 1.0), market);
        const double eu = tree(500, european(std::make_unique<PutPayoff>(100.0), 1.0), market);
        EXPECT_GE(us, eu);
        EXPECT_GE(us, std::max(100.0 - spot, 0.0));  // never below immediate exercise
    }
    // In the money with r > 0 the early exercise premium is strictly positive.
    const MarketData itm(80.0, 0.06, 0.0, 0.3);
    EXPECT_GT(tree(500, american(std::make_unique<PutPayoff>(100.0), 1.0), itm) -
                  tree(500, european(std::make_unique<PutPayoff>(100.0), 1.0), itm),
              0.1);
}

TEST(Binomial, AmericanCallWithoutDividendsEqualsEuropean) {
    // Merton: with q = 0 and r >= 0, C >= S - K e^{-rT} > S - K, so early exercise is never
    // optimal and the max() never binds: the two trees give the same number.
    for (double spot : {80.0, 100.0, 140.0}) {
        const MarketData market(spot, 0.05, 0.0, 0.25);
        const double us = tree(500, american(std::make_unique<CallPayoff>(100.0), 1.0), market);
        const double eu = tree(500, european(std::make_unique<CallPayoff>(100.0), 1.0), market);
        EXPECT_NEAR(us, eu, 1e-12 * eu);
    }
}

TEST(Binomial, AmericanCallWithDividendsCanExceedEuropean) {
    // A large dividend yield makes holding the stock attractive: exercise early to get it.
    const MarketData market(150.0, 0.02, 0.10, 0.2);
    const double us = tree(500, american(std::make_unique<CallPayoff>(100.0), 2.0), market);
    const double eu = tree(500, european(std::make_unique<CallPayoff>(100.0), 2.0), market);
    EXPECT_GT(us - eu, 1.0);
    EXPECT_NEAR(us, 50.0, 1e-9);  // deep ITM: exercise immediately
}

TEST(Binomial, DeepInTheMoneyAmericanPutIsExercisedImmediately) {
    const MarketData market(20.0, 0.10, 0.0, 0.2);
    EXPECT_NEAR(tree(500, american(std::make_unique<PutPayoff>(100.0), 1.0), market), 80.0, 1e-9);
}

TEST(Binomial, AmericanPutIncreasesWithMaturity) {
    // More time = more exercise opportunities (not true for European puts in general).
    const MarketData market(100.0, 0.08, 0.0, 0.25);
    double previous = 0.0;
    for (double t : {0.25, 0.5, 1.0, 2.0, 5.0}) {
        const double value = tree(400, american(std::make_unique<PutPayoff>(100.0), t), market);
        EXPECT_GT(value, previous);
        previous = value;
    }
}

// ---------------------------------------------------------------------------
// Edge cases and validation
// ---------------------------------------------------------------------------
TEST(Binomial, ZeroMaturityIsIntrinsic) {
    const MarketData market(110.0, 0.05, 0.0, 0.2);
    EXPECT_DOUBLE_EQ(tree(10, american(std::make_unique<CallPayoff>(100.0), 0.0), market), 10.0);
}

TEST(Binomial, OneStepTreeByHand) {
    // S=100, u=e^{0.2}, d=e^{-0.2}, r=q=0: p = (1-d)/(u-d); call K=100 pays 100(u-1) up.
    const MarketData market(100.0, 0.0, 0.0, 0.2);
    const double u = std::exp(0.2), d = std::exp(-0.2);
    const double p = (1.0 - d) / (u - d);
    EXPECT_NEAR(tree(1, european(std::make_unique<CallPayoff>(100.0), 1.0), market),
                p * 100.0 * (u - 1.0), 1e-12);
}

TEST(Binomial, RejectsInvalidSetups) {
    EXPECT_THROW(BinomialPricer(0), std::invalid_argument);
    const Option call = european(std::make_unique<CallPayoff>(100.0), 1.0);
    EXPECT_THROW((void)tree(10, call, MarketData(100.0, 0.05, 0.0, 0.0)), std::invalid_argument);
    // sigma = 1%, r = 50%, one step: e^{(r-q)dt} > u, so p > 1 (arbitrage in the tree).
    EXPECT_THROW((void)tree(1, call, MarketData(100.0, 0.5, 0.0, 0.01)), std::domain_error);
    EXPECT_NO_THROW((void)tree(5000, call, MarketData(100.0, 0.5, 0.0, 0.01)));
}
