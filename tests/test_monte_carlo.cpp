#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>

#include "pricer/analytic_pricer.hpp"
#include "pricer/monte_carlo_pricer.hpp"

using namespace pricer;

namespace {

const MarketData kHullMarket(42.0, 0.10, 0.0, 0.20);

double analytic(const Option& option, const MarketData& market) {
    return AnalyticPricer().price(option, market).price;
}

bool insideCi(const PricingResult& r, double value) {
    return r.confidenceLow <= value && value <= r.confidenceHigh;
}

// E[S_T^2] payoff: not handled by the analytic pricer, but MC prices anything.
class SquaredPayoff final : public Payoff {
public:
    double operator()(double spot) const override { return spot * spot; }
    std::unique_ptr<Payoff> clone() const override { return std::make_unique<SquaredPayoff>(); }
};

}  // namespace

TEST(MonteCarlo, HullCallAndPutInsideConfidenceInterval) {
    const MonteCarloPricer mc(200000);
    const Option call(std::make_unique<CallPayoff>(40.0), 0.5);
    const Option put(std::make_unique<PutPayoff>(40.0), 0.5);

    const PricingResult c = mc.price(call, kHullMarket);
    const PricingResult p = mc.price(put, kHullMarket);
    EXPECT_TRUE(insideCi(c, analytic(call, kHullMarket)))
        << c.price << " +/- " << 1.96 * c.standardError;
    EXPECT_TRUE(insideCi(p, analytic(put, kHullMarket)))
        << p.price << " +/- " << 1.96 * p.standardError;
    EXPECT_LT(c.standardError, 0.02);
    EXPECT_NEAR(c.confidenceHigh - c.confidenceLow, 2.0 * MonteCarloPricer::kZ95 * c.standardError,
                1e-12);
}

TEST(MonteCarlo, VariousPayoffsAndDividendsInsideConfidenceInterval) {
    const MonteCarloPricer mc(200000, 7);
    const MarketData market(100.0, 0.03, 0.02, 0.35);
    const Option options[] = {
        Option(std::make_unique<CallPayoff>(90.0), 2.0),
        Option(std::make_unique<PutPayoff>(110.0), 0.25),
        Option(std::make_unique<DigitalCallPayoff>(100.0), 1.0),
        Option(std::make_unique<DigitalPutPayoff>(95.0), 1.0),
    };
    for (const Option& option : options) {
        const PricingResult r = mc.price(option, market);
        EXPECT_TRUE(insideCi(r, analytic(option, market)))
            << r.price << " vs " << analytic(option, market);
    }
}

TEST(MonteCarlo, ConfidenceIntervalCoverageIsAbout95Percent) {
    // Over 200 independent seeds the true price should fall inside the 95% CI about
    // 190 times. Binomial(200, 0.95) has standard deviation ~3.1, so [178, 199] is a
    // +/- 4 sigma band: this checks that the SE is neither under- nor over-estimated.
    const Option call(std::make_unique<CallPayoff>(40.0), 0.5);
    const double truth = analytic(call, kHullMarket);
    int covered = 0;
    for (std::uint64_t seed = 1; seed <= 200; ++seed) {
        covered += insideCi(MonteCarloPricer(5000, seed).price(call, kHullMarket), truth) ? 1 : 0;
    }
    EXPECT_GE(covered, 178);
    EXPECT_LE(covered, 199);
}

TEST(MonteCarlo, StandardErrorDecreasesAsOneOverSqrtN) {
    const Option call(std::make_unique<CallPayoff>(40.0), 0.5);
    const double se1 = MonteCarloPricer(25000, 3).price(call, kHullMarket).standardError;
    const double se4 = MonteCarloPricer(100000, 3).price(call, kHullMarket).standardError;
    const double se16 = MonteCarloPricer(400000, 3).price(call, kHullMarket).standardError;
    EXPECT_NEAR(se4 / se1, 0.5, 0.02);
    EXPECT_NEAR(se16 / se4, 0.5, 0.02);
}

TEST(MonteCarlo, ReproducibleWithSameSeedDifferentOtherwise) {
    const Option call(std::make_unique<CallPayoff>(40.0), 0.5);
    const MonteCarloPricer a(10000, 99), b(10000, 99), c(10000, 100);
    EXPECT_EQ(a.price(call, kHullMarket).price, b.price(call, kHullMarket).price);
    EXPECT_EQ(a.price(call, kHullMarket).price, a.price(call, kHullMarket).price);  // const call
    EXPECT_NE(a.price(call, kHullMarket).price, c.price(call, kHullMarket).price);
}

TEST(MonteCarlo, GenericPayoffSecondMoment) {
    // e^{-rT} E[S_T^2] = S^2 exp((r - 2q + sigma^2) T)
    const MarketData market(100.0, 0.05, 0.01, 0.2);
    const Option option(std::make_unique<SquaredPayoff>(), 1.0);
    const double exact = 100.0 * 100.0 * std::exp((0.05 - 2.0 * 0.01 + 0.04) * 1.0);
    const PricingResult r = MonteCarloPricer(200000).price(option, market);
    EXPECT_TRUE(insideCi(r, exact)) << r.price << " vs " << exact;
}

TEST(MonteCarlo, DeterministicCasesHaveZeroError) {
    const Option call(std::make_unique<CallPayoff>(100.0), 1.0);
    const MarketData noVol(110.0, 0.05, 0.0, 0.0);
    const PricingResult r = MonteCarloPricer(1000).price(call, noVol);
    EXPECT_NEAR(r.price, analytic(call, noVol), 1e-12);
    EXPECT_EQ(r.standardError, 0.0);

    const Option expired(std::make_unique<CallPayoff>(100.0), 0.0);
    const PricingResult e = MonteCarloPricer(1000).price(expired, kHullMarket.withSpot(110.0));
    EXPECT_DOUBLE_EQ(e.price, 10.0);
    EXPECT_EQ(e.standardError, 0.0);
}

TEST(MonteCarlo, RejectsInvalidSetup) {
    EXPECT_THROW(MonteCarloPricer(1), std::invalid_argument);
    const Option american(std::make_unique<PutPayoff>(40.0), 0.5, ExerciseStyle::American);
    EXPECT_THROW((void)MonteCarloPricer(1000).price(american, kHullMarket), std::invalid_argument);
}
