#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "pricer/black_scholes.hpp"
#include "pricer/exotic_formulas.hpp"
#include "pricer/path_dependent_pricer.hpp"

using namespace pricer;

namespace {

const MarketData kMarket(100.0, 0.05, 0.02, 0.25);
const black_scholes::Inputs kInputs{100.0, 100.0, 0.05, 0.02, 0.25, 1.0};
constexpr double kMaturity = 1.0;
constexpr std::size_t kSteps = 50;
constexpr std::size_t kPaths = 100000;

bool insideCi(const PricingResult& r, double value) {
    return r.confidenceLow <= value && value <= r.confidenceHigh;
}

}  // namespace

// ---------------------------------------------------------------------------
// Payoffs on hand-made paths
// ---------------------------------------------------------------------------
TEST(PathPayoffs, ValuesOnAKnownPath) {
    const std::vector<double> path = {90.0, 110.0, 120.0, 100.0};  // mean 105, last 100
    EXPECT_DOUBLE_EQ(ArithmeticAsianPayoff(OptionType::Call, 100.0)(path), 5.0);
    EXPECT_DOUBLE_EQ(ArithmeticAsianPayoff(OptionType::Put, 110.0)(path), 5.0);
    const double geometric = std::pow(90.0 * 110.0 * 120.0 * 100.0, 0.25);
    EXPECT_NEAR(GeometricAsianPayoff(OptionType::Call, 100.0)(path), geometric - 100.0, 1e-12);
    EXPECT_DOUBLE_EQ(TerminalPayoff(std::make_unique<PutPayoff>(105.0))(path), 5.0);

    // Max 120, min 90.
    auto put = [&](BarrierType type, double level) {
        return BarrierPayoff(OptionType::Put, 105.0, type, level)(path);
    };
    EXPECT_DOUBLE_EQ(put(BarrierType::UpAndOut, 120.0), 0.0);    // touched at 120
    EXPECT_DOUBLE_EQ(put(BarrierType::UpAndOut, 121.0), 5.0);    // never touched
    EXPECT_DOUBLE_EQ(put(BarrierType::UpAndIn, 120.0), 5.0);
    EXPECT_DOUBLE_EQ(put(BarrierType::DownAndOut, 90.0), 0.0);   // touched at 90
    EXPECT_DOUBLE_EQ(put(BarrierType::DownAndIn, 89.0), 0.0);    // never touched
}

TEST(PathPayoffs, GeometricMeanDoesNotOverflow) {
    // 400 fixings around 1e3: the naive product would be 1e1200 = inf.
    const std::vector<double> path(400, 1000.0);
    EXPECT_NEAR(GeometricAsianPayoff(OptionType::Call, 900.0)(path), 100.0, 1e-9);
}

TEST(PathPayoffs, RejectInvalidParameters) {
    EXPECT_THROW(ArithmeticAsianPayoff(OptionType::Call, 0.0), std::invalid_argument);
    EXPECT_THROW(BarrierPayoff(OptionType::Call, 100.0, BarrierType::UpAndOut, -1.0),
                 std::invalid_argument);
    EXPECT_THROW(TerminalPayoff(nullptr), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Path simulation sanity
// ---------------------------------------------------------------------------
TEST(PathDependentPricer, TerminalPayoffMatchesBlackScholes) {
    const PathDependentPricer pricer(kPaths, kSteps, 1);
    const TerminalPayoff call(std::make_unique<CallPayoff>(100.0));
    const TerminalPayoff put(std::make_unique<PutPayoff>(90.0));
    EXPECT_TRUE(insideCi(pricer.price(call, kMaturity, kMarket),
                         black_scholes::vanillaPrice(OptionType::Call, kInputs)));
    black_scholes::Inputs putInputs = kInputs;
    putInputs.strike = 90.0;
    EXPECT_TRUE(insideCi(pricer.price(put, kMaturity, kMarket),
                         black_scholes::vanillaPrice(OptionType::Put, putInputs)));
}

TEST(PathDependentPricer, RejectsInvalidSetup) {
    EXPECT_THROW(PathDependentPricer(1, 10), std::invalid_argument);
    EXPECT_THROW(PathDependentPricer(100, 0), std::invalid_argument);
    const TerminalPayoff call(std::make_unique<CallPayoff>(100.0));
    EXPECT_THROW((void)PathDependentPricer(100, 10).price(call, 0.0, kMarket),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Asian options
// ---------------------------------------------------------------------------
TEST(Asian, GeometricClosedFormReducesToBlackScholesWithOneFixing) {
    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        EXPECT_NEAR(closed_form::geometricAsianPrice(type, kInputs, 1),
                    black_scholes::vanillaPrice(type, kInputs), 1e-12);
    }
}

TEST(Asian, GeometricMonteCarloMatchesClosedForm) {
    const PathDependentPricer pricer(kPaths, kSteps, 2);
    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        const PricingResult mc =
            pricer.price(GeometricAsianPayoff(type, 100.0), kMaturity, kMarket);
        EXPECT_TRUE(insideCi(mc, closed_form::geometricAsianPrice(type, kInputs, kSteps)));
    }
}

TEST(Asian, ArithmeticAboveGeometricBelowVanilla) {
    // AM >= GM path by path, so the arithmetic call is worth more than the geometric one.
    // Averaging reduces volatility, so the Asian call is cheaper than the vanilla call.
    const PathDependentPricer pricer(kPaths, kSteps, 3);
    const PricingResult arithmetic =
        pricer.price(ArithmeticAsianPayoff(OptionType::Call, 100.0), kMaturity, kMarket);
    EXPECT_GT(arithmetic.confidenceLow,
              closed_form::geometricAsianPrice(OptionType::Call, kInputs, kSteps));
    EXPECT_LT(arithmetic.confidenceHigh, black_scholes::vanillaPrice(OptionType::Call, kInputs));
}

TEST(Asian, GeometricControlVariateIsUnbiasedAndCutsVarianceByOver100) {
    const PathDependentPricer pricer(kPaths, kSteps, 4);
    const ArithmeticAsianPayoff target(OptionType::Call, 100.0);
    const GeometricAsianPayoff control(OptionType::Call, 100.0);
    const double controlPrice = closed_form::geometricAsianPrice(OptionType::Call, kInputs, kSteps);

    const PricingResult plain = pricer.price(target, kMaturity, kMarket);
    const PricingResult withControl =
        pricer.priceWithControlVariate(target, control, controlPrice, kMaturity, kMarket);
    EXPECT_TRUE(insideCi(plain, withControl.price));  // consistent estimates
    const double factor = std::pow(plain.standardError / withControl.standardError, 2.0);
    EXPECT_GT(factor, 100.0);

    // A perfect control (the target itself) removes all the variance.
    const PricingResult self = pricer.priceWithControlVariate(
        control, control, controlPrice, kMaturity, kMarket);
    EXPECT_NEAR(self.price, controlPrice, 1e-10);
    EXPECT_LT(self.standardError, 1e-10);
}

// ---------------------------------------------------------------------------
// Barrier options
// ---------------------------------------------------------------------------
TEST(Barrier, InPlusOutEqualsVanillaPathByPath) {
    // Same seed, so the same paths: the identity holds to rounding, not just statistically.
    const PathDependentPricer pricer(20000, kSteps, 5);
    const TerminalPayoff vanillaCall(std::make_unique<CallPayoff>(100.0));
    const double vanilla = pricer.price(vanillaCall, kMaturity, kMarket).price;
    const struct {
        BarrierType out, in;
        double level;
    } pairs[] = {{BarrierType::DownAndOut, BarrierType::DownAndIn, 90.0},
                 {BarrierType::UpAndOut, BarrierType::UpAndIn, 120.0}};
    for (const auto& p : pairs) {
        const double out =
            pricer.price(BarrierPayoff(OptionType::Call, 100.0, p.out, p.level), kMaturity, kMarket)
                .price;
        const double in =
            pricer.price(BarrierPayoff(OptionType::Call, 100.0, p.in, p.level), kMaturity, kMarket)
                .price;
        EXPECT_NEAR(out + in, vanilla, 1e-10);
    }
}

TEST(Barrier, DiscreteMonitoringIsWorthMoreThanContinuousForKnockOut) {
    // Crossings between monitoring dates are missed, so a discretely monitored knock-out
    // survives more often than the continuous one.
    const double barrier = 90.0;
    const double continuous = closed_form::downAndOutCallPrice(kInputs, barrier);
    const PricingResult discrete = PathDependentPricer(kPaths, kSteps, 6).price(
        BarrierPayoff(OptionType::Call, 100.0, BarrierType::DownAndOut, barrier), kMaturity,
        kMarket);
    EXPECT_GT(discrete.confidenceLow, continuous);
}

TEST(Barrier, BroadieGlassermanKouCorrectionMatchesDiscreteMonitoring) {
    // Continuous formula evaluated at the shifted barrier H exp(-0.5826 sigma sqrt(dt))
    // reproduces the discretely monitored price.
    const double barrier = 90.0;
    for (std::size_t steps : {12U, 50U}) {
        const double dt = kMaturity / static_cast<double>(steps);
        const double corrected = closed_form::downAndOutCallPrice(
            kInputs, closed_form::continuityCorrectedBarrier(barrier, 0.25, dt, false));
        const PricingResult discrete = PathDependentPricer(kPaths, steps, 7).price(
            BarrierPayoff(OptionType::Call, 100.0, BarrierType::DownAndOut, barrier), kMaturity,
            kMarket);
        EXPECT_NEAR(discrete.price, corrected, 3.0 * discrete.standardError) << "M=" << steps;
    }
}

TEST(Barrier, ClosedFormLimits) {
    // Barrier far below the spot: knock-out ~ vanilla. In + out = vanilla by construction.
    EXPECT_NEAR(closed_form::downAndOutCallPrice(kInputs, 1.0),
                black_scholes::vanillaPrice(OptionType::Call, kInputs), 1e-10);
    EXPECT_GT(closed_form::downAndInCallPrice(kInputs, 90.0), 0.0);
    EXPECT_THROW((void)closed_form::downAndInCallPrice(kInputs, 110.0), std::invalid_argument);
}
