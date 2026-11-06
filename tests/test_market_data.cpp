#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>

#include "pricer/market_data.hpp"

using namespace pricer;

TEST(MarketData, StoresInputs) {
    const MarketData m(42.0, 0.10, 0.02, 0.20);
    EXPECT_DOUBLE_EQ(m.spot(), 42.0);
    EXPECT_DOUBLE_EQ(m.rate(), 0.10);
    EXPECT_DOUBLE_EQ(m.dividend(), 0.02);
    EXPECT_DOUBLE_EQ(m.volatility(), 0.20);
}

TEST(MarketData, AcceptsNegativeRatesAndZeroVolatility) {
    EXPECT_NO_THROW(MarketData(100.0, -0.005, -0.01, 0.0));
}

TEST(MarketData, RejectsInvalidInputs) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(MarketData(0.0, 0.05, 0.0, 0.2), std::invalid_argument);
    EXPECT_THROW(MarketData(-1.0, 0.05, 0.0, 0.2), std::invalid_argument);
    EXPECT_THROW(MarketData(100.0, 0.05, 0.0, -0.2), std::invalid_argument);
    EXPECT_THROW(MarketData(100.0, nan, 0.0, 0.2), std::invalid_argument);
}

TEST(MarketData, WithBumpsReturnModifiedCopy) {
    const MarketData m(100.0, 0.05, 0.01, 0.2);
    const MarketData up = m.withSpot(101.0).withVolatility(0.25).withRate(0.06);
    EXPECT_DOUBLE_EQ(up.spot(), 101.0);
    EXPECT_DOUBLE_EQ(up.volatility(), 0.25);
    EXPECT_DOUBLE_EQ(up.rate(), 0.06);
    EXPECT_DOUBLE_EQ(up.dividend(), 0.01);
    EXPECT_DOUBLE_EQ(m.spot(), 100.0);  // original untouched
    EXPECT_THROW((void)m.withSpot(-5.0), std::invalid_argument);
}
