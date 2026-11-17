#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>

#include "pricer/analytic_pricer.hpp"
#include "pricer/monte_carlo_pricer.hpp"
#include "pricer/statistics.hpp"

using namespace pricer;

namespace {

const MarketData kHullMarket(42.0, 0.10, 0.0, 0.20);
constexpr std::size_t kPaths = 100000;

double analytic(const Option& option, const MarketData& market) {
    return AnalyticPricer().price(option, market).price;
}

bool insideCi(const PricingResult& r, double value) {
    return r.confidenceLow <= value && value <= r.confidenceHigh;
}

/// Variance reduction factor at a fixed number of payoff evaluations (i.e. at equal cost):
/// Var(plain) / Var(method) = (SE_plain / SE_method)^2.
double reductionFactor(const Option& option, const MarketData& market, VarianceReduction vr,
                       std::uint64_t seed = 11) {
    const double plain = MonteCarloPricer(kPaths, seed).price(option, market).standardError;
    const double reduced = MonteCarloPricer(kPaths, seed, vr).price(option, market).standardError;
    return (plain * plain) / (reduced * reduced);
}

// Long straddle |S_T - K|: a symmetric (non-monotone) payoff in Z near the forward.
class StraddlePayoff final : public Payoff {
public:
    explicit StraddlePayoff(double strike) : strike_(strike) {}
    double operator()(double spot) const override { return std::abs(spot - strike_); }
    std::unique_ptr<Payoff> clone() const override {
        return std::make_unique<StraddlePayoff>(*this);
    }

private:
    double strike_;
};

}  // namespace

// ---------------------------------------------------------------------------
// RunningCovariance
// ---------------------------------------------------------------------------
TEST(RunningCovariance, KnownValues) {
    RunningCovariance c;
    // y = 2x + 1 exactly: Cov = 2 Var(x), Var(y) = 4 Var(x).
    for (double x : {1.0, 2.0, 4.0, 7.0}) c.add(x, 2.0 * x + 1.0);
    EXPECT_DOUBLE_EQ(c.meanX(), 3.5);
    EXPECT_DOUBLE_EQ(c.meanY(), 8.0);
    EXPECT_DOUBLE_EQ(c.varianceX(), 7.0);
    EXPECT_DOUBLE_EQ(c.covariance(), 14.0);
    EXPECT_DOUBLE_EQ(c.varianceY(), 28.0);
}

// ---------------------------------------------------------------------------
// Unbiasedness: the Black-Scholes price lies inside the 95% CI for every method
// ---------------------------------------------------------------------------
class VarianceReductionMethods : public ::testing::TestWithParam<VarianceReduction> {};

TEST_P(VarianceReductionMethods, BlackScholesInsideConfidenceInterval) {
    const MarketData market(100.0, 0.03, 0.02, 0.35);
    const Option options[] = {
        Option(std::make_unique<CallPayoff>(40.0), 0.5),
        Option(std::make_unique<CallPayoff>(90.0), 2.0),
        Option(std::make_unique<PutPayoff>(110.0), 0.25),
        Option(std::make_unique<DigitalCallPayoff>(100.0), 1.0),
    };
    const MonteCarloPricer mc(kPaths, 5, GetParam());
    EXPECT_TRUE(insideCi(mc.price(options[0], kHullMarket), analytic(options[0], kHullMarket)));
    for (int i = 1; i < 4; ++i) {
        const PricingResult r = mc.price(options[i], market);
        EXPECT_TRUE(insideCi(r, analytic(options[i], market)))
            << "option " << i << ": " << r.price << " vs " << analytic(options[i], market);
    }
}

TEST_P(VarianceReductionMethods, CoverageIsAbout95Percent) {
    // Checks both that the estimator is (asymptotically) unbiased and that its SE is right.
    // Includes the control variate, whose beta is estimated on the same paths (O(1/N) bias).
    const Option call(std::make_unique<CallPayoff>(40.0), 0.5);
    const double truth = analytic(call, kHullMarket);
    int covered = 0;
    for (std::uint64_t seed = 1; seed <= 200; ++seed) {
        covered += insideCi(MonteCarloPricer(2000, seed, GetParam()).price(call, kHullMarket),
                            truth) ? 1 : 0;
    }
    EXPECT_GE(covered, 178);
    EXPECT_LE(covered, 199);
}

TEST_P(VarianceReductionMethods, ReproducibleAndZeroErrorWhenDeterministic) {
    const Option call(std::make_unique<CallPayoff>(100.0), 1.0);
    const MonteCarloPricer mc(1000, 3, GetParam());
    EXPECT_EQ(mc.price(call, kHullMarket.withSpot(100.0)).price,
              mc.price(call, kHullMarket.withSpot(100.0)).price);

    const MarketData noVol(110.0, 0.05, 0.0, 0.0);
    const PricingResult r = mc.price(call, noVol);
    EXPECT_NEAR(r.price, analytic(call, noVol), 1e-12);
    EXPECT_EQ(r.standardError, 0.0);
}

INSTANTIATE_TEST_SUITE_P(All, VarianceReductionMethods,
                         ::testing::Values(VarianceReduction::None, VarianceReduction::Antithetic,
                                           VarianceReduction::ControlVariate));

// ---------------------------------------------------------------------------
// Variance reduction factors (at equal number of payoff evaluations)
// ---------------------------------------------------------------------------
TEST(VarianceReduction, AntitheticHelpsForMonotonePayoffs) {
    const Option call(std::make_unique<CallPayoff>(40.0), 0.5);
    const Option put(std::make_unique<PutPayoff>(40.0), 0.5);
    EXPECT_GT(reductionFactor(call, kHullMarket, VarianceReduction::Antithetic), 2.0);
    EXPECT_GT(reductionFactor(put, kHullMarket, VarianceReduction::Antithetic), 1.2);
}

TEST(VarianceReduction, AntitheticHurtsForSymmetricPayoff) {
    // For a straddle struck at the forward, h(Z) and h(-Z) are positively correlated:
    // pairing them wastes half the budget and the variance goes UP.
    const MarketData market(100.0, 0.0, 0.0, 0.2);
    const Option straddle(std::make_unique<StraddlePayoff>(100.0), 1.0);
    EXPECT_LT(reductionFactor(straddle, market, VarianceReduction::Antithetic), 1.0);
}

TEST(VarianceReduction, ControlVariateGainGrowsWithMoneyness) {
    // Var(h - beta X) = Var(h) (1 - rho^2): the gain is 1 / (1 - rho^2), large when the
    // call behaves like the underlying (deep ITM), small when it rarely pays (deep OTM).
    const MarketData market(100.0, 0.05, 0.0, 0.2);
    const Option deepItm(std::make_unique<CallPayoff>(60.0), 1.0);
    const Option atm(std::make_unique<CallPayoff>(100.0), 1.0);
    const Option deepOtm(std::make_unique<CallPayoff>(160.0), 1.0);
    const double itmGain = reductionFactor(deepItm, market, VarianceReduction::ControlVariate);
    const double atmGain = reductionFactor(atm, market, VarianceReduction::ControlVariate);
    const double otmGain = reductionFactor(deepOtm, market, VarianceReduction::ControlVariate);
    EXPECT_GT(itmGain, 50.0);
    EXPECT_GT(atmGain, 2.5);
    EXPECT_GT(itmGain, atmGain);
    EXPECT_GT(atmGain, otmGain);
    EXPECT_GE(otmGain, 1.0);
}

TEST(VarianceReduction, RejectsTooFewPathsForAntithetic) {
    EXPECT_THROW(MonteCarloPricer(3, 1, VarianceReduction::Antithetic), std::invalid_argument);
    EXPECT_NO_THROW(MonteCarloPricer(4, 1, VarianceReduction::Antithetic));
}
