#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>
#include <utility>

#include "pricer/option.hpp"

using namespace pricer;

TEST(Option, StoresContractTerms) {
    const Option opt(std::make_unique<CallPayoff>(40.0), 0.5, ExerciseStyle::American);
    EXPECT_DOUBLE_EQ(opt.maturity(), 0.5);
    EXPECT_EQ(opt.exercise(), ExerciseStyle::American);
    EXPECT_DOUBLE_EQ(opt.payoffAt(42.0), 2.0);
}

TEST(Option, DefaultsToEuropean) {
    const Option opt(std::make_unique<PutPayoff>(40.0), 1.0);
    EXPECT_EQ(opt.exercise(), ExerciseStyle::European);
}

TEST(Option, CopyIsDeep) {
    const Option a(std::make_unique<PutPayoff>(40.0), 1.0);
    const Option b(a);  // NOLINT(performance-unnecessary-copy-initialization)
    EXPECT_NE(&a.payoff(), &b.payoff());
    EXPECT_DOUBLE_EQ(b.payoffAt(30.0), 10.0);

    Option c(std::make_unique<CallPayoff>(1.0), 2.0);
    c = a;
    EXPECT_NE(&c.payoff(), &a.payoff());
    EXPECT_DOUBLE_EQ(c.payoffAt(30.0), 10.0);
    EXPECT_DOUBLE_EQ(c.maturity(), 1.0);
}

TEST(Option, MoveTransfersOwnership) {
    Option a(std::make_unique<CallPayoff>(10.0), 1.0);
    const Payoff* raw = &a.payoff();
    const Option b(std::move(a));
    EXPECT_EQ(&b.payoff(), raw);
}

TEST(Option, RejectsInvalidInputs) {
    EXPECT_THROW(Option(nullptr, 1.0), std::invalid_argument);
    EXPECT_THROW(Option(std::make_unique<CallPayoff>(10.0), -0.1), std::invalid_argument);
}

TEST(Option, AcceptsZeroMaturity) {
    EXPECT_NO_THROW(Option(std::make_unique<CallPayoff>(10.0), 0.0));
}
