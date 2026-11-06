#include <gtest/gtest.h>

#include <limits>
#include <memory>
#include <stdexcept>

#include "pricer/payoff.hpp"

using namespace pricer;

TEST(Payoff, CallIsMaxOfSpotMinusStrikeAndZero) {
    const CallPayoff call(100.0);
    EXPECT_DOUBLE_EQ(call(120.0), 20.0);
    EXPECT_DOUBLE_EQ(call(100.0), 0.0);
    EXPECT_DOUBLE_EQ(call(80.0), 0.0);
}

TEST(Payoff, PutIsMaxOfStrikeMinusSpotAndZero) {
    const PutPayoff put(100.0);
    EXPECT_DOUBLE_EQ(put(80.0), 20.0);
    EXPECT_DOUBLE_EQ(put(100.0), 0.0);
    EXPECT_DOUBLE_EQ(put(120.0), 0.0);
}

TEST(Payoff, CallMinusPutIsForwardPayoff) {
    // Model-free identity behind put-call parity: (S-K)+ - (K-S)+ = S - K.
    const CallPayoff call(100.0);
    const PutPayoff put(100.0);
    for (double s : {50.0, 99.0, 100.0, 101.0, 250.0}) {
        EXPECT_DOUBLE_EQ(call(s) - put(s), s - 100.0);
    }
}

TEST(Payoff, DigitalsPayOneOnTheirSideOfTheStrike) {
    const DigitalCallPayoff dcall(100.0);
    const DigitalPutPayoff dput(100.0);
    EXPECT_DOUBLE_EQ(dcall(101.0), 1.0);
    EXPECT_DOUBLE_EQ(dcall(99.0), 0.0);
    EXPECT_DOUBLE_EQ(dput(99.0), 1.0);
    EXPECT_DOUBLE_EQ(dput(101.0), 0.0);
    // At S = K exactly both pay 0 (measure-zero event under a continuous law).
    EXPECT_DOUBLE_EQ(dcall(100.0), 0.0);
    EXPECT_DOUBLE_EQ(dput(100.0), 0.0);
}

TEST(Payoff, CloneKeepsDynamicTypeAndValue) {
    const std::unique_ptr<Payoff> original = std::make_unique<PutPayoff>(50.0);
    const std::unique_ptr<Payoff> copy = original->clone();
    ASSERT_NE(copy, nullptr);
    EXPECT_NE(copy.get(), original.get());
    EXPECT_NE(dynamic_cast<const PutPayoff*>(copy.get()), nullptr);
    EXPECT_DOUBLE_EQ((*copy)(40.0), 10.0);
}

TEST(Payoff, RejectsInvalidStrike) {
    EXPECT_THROW(CallPayoff(0.0), std::invalid_argument);
    EXPECT_THROW(PutPayoff(-1.0), std::invalid_argument);
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(DigitalCallPayoff{inf}, std::invalid_argument);
}
