#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "pricer/analytic_pricer.hpp"
#include "pricer/black_scholes.hpp"

using namespace pricer;

namespace {

const MarketData kHullMarket(42.0, 0.10, 0.0, 0.20);

// A payoff the analytic pricer does not know about.
class SquaredPayoff final : public Payoff {
public:
    double operator()(double spot) const override { return spot * spot; }
    std::unique_ptr<Payoff> clone() const override { return std::make_unique<SquaredPayoff>(); }
};

}  // namespace

TEST(AnalyticPricer, PricesHullCallAndPutThroughOptionInterface) {
    const AnalyticPricer pricer;
    const Option call(std::make_unique<CallPayoff>(40.0), 0.5);
    const Option put(std::make_unique<PutPayoff>(40.0), 0.5);

    const PricingResult c = pricer.price(call, kHullMarket);
    const PricingResult p = pricer.price(put, kHullMarket);
    EXPECT_NEAR(c.price, 4.76, 0.005);
    EXPECT_NEAR(p.price, 0.81, 0.005);
    EXPECT_EQ(c.standardError, 0.0);
    EXPECT_EQ(c.confidenceLow, c.price);
    EXPECT_EQ(c.confidenceHigh, c.price);
}

TEST(AnalyticPricer, WorksThroughBaseClassReference) {
    const AnalyticPricer analytic;
    const Pricer& pricer = analytic;
    const Option digital(std::make_unique<DigitalCallPayoff>(40.0), 0.5);
    const double expected = black_scholes::digitalPrice(
        OptionType::Call, {42.0, 40.0, 0.10, 0.0, 0.20, 0.5});
    EXPECT_DOUBLE_EQ(pricer.price(digital, kHullMarket).price, expected);
}

TEST(AnalyticPricer, RejectsAmericanExercise) {
    const AnalyticPricer pricer;
    const Option american(std::make_unique<PutPayoff>(40.0), 0.5, ExerciseStyle::American);
    EXPECT_THROW((void)pricer.price(american, kHullMarket), std::invalid_argument);
}

TEST(AnalyticPricer, RejectsUnknownPayoff) {
    const AnalyticPricer pricer;
    const Option exotic(std::make_unique<SquaredPayoff>(), 0.5);
    EXPECT_THROW((void)pricer.price(exotic, kHullMarket), std::invalid_argument);
}
