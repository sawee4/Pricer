#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>

#include "pricer/analytic_pricer.hpp"
#include "pricer/monte_carlo_greeks.hpp"
#include "pricer/statistics.hpp"

using namespace pricer;
using namespace pricer::mc_greeks;

namespace {

const MarketData kMarket(100.0, 0.05, 0.02, 0.25);
constexpr std::size_t kPaths = 200000;

Greeks analyticGreeks(const Option& option, const MarketData& market) {
    return AnalyticPricer().greeks(option, market);
}

void expectRelNear(double actual, double expected, double rel) {
    EXPECT_NEAR(actual, expected, rel * std::abs(expected)) << "expected " << expected;
}

}  // namespace

// ---------------------------------------------------------------------------
// Bump-and-revalue with common random numbers
// ---------------------------------------------------------------------------
TEST(BumpAndRevalue, CommonRandomNumbersMatchAnalyticGreeks) {
    const Option call(std::make_unique<CallPayoff>(105.0), 1.0);
    const Option put(std::make_unique<PutPayoff>(95.0), 0.5);
    for (const Option* option : {&call, &put}) {
        const Greeks exact = analyticGreeks(*option, kMarket);
        const Greeks mc = bumpAndRevalue(MonteCarloPricer(kPaths, 1), *option, kMarket);
        EXPECT_NEAR(mc.delta, exact.delta, 0.005);
        expectRelNear(mc.gamma, exact.gamma, 0.05);
        expectRelNear(mc.vega, exact.vega, 0.02);
        expectRelNear(mc.theta, exact.theta, 0.03);
        expectRelNear(mc.rho, exact.rho, 0.02);
    }
}

TEST(BumpAndRevalue, CommonRandomNumbersCutDeltaNoise) {
    // Standard deviation of the delta estimator across 40 seeds, with and without CRN.
    // Independent draws: Var ~ 2 Var(V) / (4 h^2 N), huge for h = 1% of S.
    // CRN: the noise of V(S+h) and V(S-h) cancels; Var ~ Var(pathwise delta) / N.
    const Option call(std::make_unique<CallPayoff>(105.0), 1.0);
    BumpSettings independent;
    independent.commonRandomNumbers = false;
    RunningStats crn, ind;
    for (std::uint64_t seed = 1; seed <= 40; ++seed) {
        const MonteCarloPricer pricer(20000, seed);
        crn.add(bumpAndRevalue(pricer, call, kMarket).delta);
        ind.add(bumpAndRevalue(pricer, call, kMarket, independent).delta);
    }
    EXPECT_GT(std::sqrt(ind.variance()), 10.0 * std::sqrt(crn.variance()));
}

TEST(BumpAndRevalue, SmallBumpWithCrnConvergesToPathwiseDelta) {
    // With the same draws, [h(S_T(S+h)) - h(S_T(S-h))] / 2h -> h'(S_T) dS_T/dS_0 path by path
    // as h -> 0 (except on the measure-zero kink): bump-and-revalue with CRN *is* the
    // pathwise estimator in the limit.
    const Option call(std::make_unique<CallPayoff>(105.0), 1.0);
    BumpSettings tiny;
    tiny.relativeSpotBump = 1e-7;
    const double bumped = bumpAndRevalue(MonteCarloPricer(50000, 9), call, kMarket, tiny).delta;
    EXPECT_NEAR(bumped, pathwiseDelta(call, kMarket, 50000, 9).value, 1e-6);
}

TEST(BumpAndRevalue, RejectsInvalidSettings) {
    const Option call(std::make_unique<CallPayoff>(100.0), 1.0);
    const MonteCarloPricer pricer(1000);
    BumpSettings bad;
    bad.relativeSpotBump = 0.0;
    EXPECT_THROW((void)bumpAndRevalue(pricer, call, kMarket, bad), std::invalid_argument);
    const Option shortDated(std::make_unique<CallPayoff>(100.0), 0.5 / 365.0);
    EXPECT_THROW((void)bumpAndRevalue(pricer, shortDated, kMarket), std::invalid_argument);
    EXPECT_THROW((void)bumpAndRevalue(pricer, call, kMarket.withVolatility(0.0005)),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Pathwise estimators
// ---------------------------------------------------------------------------
TEST(Pathwise, DeltaAndVegaWithinThreeStandardErrors) {
    const MarketData markets[] = {kMarket, MarketData(42.0, 0.10, 0.0, 0.20)};
    for (const MarketData& market : markets) {
        for (double strike : {0.8 * market.spot(), market.spot(), 1.2 * market.spot()}) {
            const Option call(std::make_unique<CallPayoff>(strike), 0.75);
            const Option put(std::make_unique<PutPayoff>(strike), 0.75);
            for (const Option* option : {&call, &put}) {
                const Greeks exact = analyticGreeks(*option, market);
                const Estimate delta = pathwiseDelta(*option, market, kPaths, 3);
                const Estimate vega = pathwiseVega(*option, market, kPaths, 3);
                EXPECT_NEAR(delta.value, exact.delta, 3.0 * delta.standardError);
                EXPECT_NEAR(vega.value, exact.vega, 3.0 * vega.standardError);
            }
        }
    }
}

TEST(Pathwise, LowerVarianceThanLikelihoodRatio) {
    const Option call(std::make_unique<CallPayoff>(105.0), 1.0);
    const Estimate pw = pathwiseDelta(call, kMarket, kPaths, 4);
    const Estimate lr = likelihoodRatioDelta(call, kMarket, kPaths, 4);
    EXPECT_LT(pw.standardError, lr.standardError);
}

TEST(Pathwise, RejectsNonLipschitzPayoffs) {
    // For a digital, h'(S_T) = 0 almost everywhere: the pathwise estimator would return 0.
    const Option digital(std::make_unique<DigitalCallPayoff>(100.0), 1.0);
    EXPECT_THROW((void)pathwiseDelta(digital, kMarket, 1000, 1), std::invalid_argument);
    EXPECT_THROW((void)pathwiseVega(digital, kMarket, 1000, 1), std::invalid_argument);
    const Option american(std::make_unique<PutPayoff>(100.0), 1.0, ExerciseStyle::American);
    EXPECT_THROW((void)pathwiseDelta(american, kMarket, 1000, 1), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Likelihood ratio
// ---------------------------------------------------------------------------
TEST(LikelihoodRatio, DeltaOfVanillaAndDigitalWithinThreeStandardErrors) {
    const Option call(std::make_unique<CallPayoff>(105.0), 1.0);
    const Estimate lrCall = likelihoodRatioDelta(call, kMarket, kPaths, 5);
    EXPECT_NEAR(lrCall.value, analyticGreeks(call, kMarket).delta, 3.0 * lrCall.standardError);

    // Exact digital delta by central difference of the closed-form digital price.
    const Option digital(std::make_unique<DigitalCallPayoff>(105.0), 1.0);
    const AnalyticPricer analytic;
    const double h = 1e-4;
    const double exact = (analytic.price(digital, kMarket.withSpot(100.0 + h)).price -
                          analytic.price(digital, kMarket.withSpot(100.0 - h)).price) /
                         (2.0 * h);
    const Estimate lrDigital = likelihoodRatioDelta(digital, kMarket, kPaths, 5);
    EXPECT_NEAR(lrDigital.value, exact, 3.0 * lrDigital.standardError);
}

TEST(LikelihoodRatio, RejectsDegenerateInputs) {
    const Option call(std::make_unique<CallPayoff>(100.0), 1.0);
    EXPECT_THROW((void)likelihoodRatioDelta(call, kMarket.withVolatility(0.0), 1000, 1),
                 std::invalid_argument);
    EXPECT_THROW((void)likelihoodRatioDelta(call, kMarket, 1, 1), std::invalid_argument);
}
