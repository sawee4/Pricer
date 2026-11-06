#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "pricer/normal.hpp"

using pricer::math::normCdf;
using pricer::math::normPdf;

TEST(Normal, KnownValues) {
    EXPECT_DOUBLE_EQ(normCdf(0.0), 0.5);
    EXPECT_NEAR(normCdf(1.959963984540054), 0.975, 1e-15);
    EXPECT_NEAR(normCdf(-1.0), 0.15865525393145707, 1e-15);
    EXPECT_NEAR(normPdf(0.0), 0.3989422804014327, 1e-16);
}

TEST(Normal, Symmetry) {
    for (double x : {0.1, 0.5, 1.0, 2.5, 5.0}) {
        EXPECT_NEAR(normCdf(x) + normCdf(-x), 1.0, 1e-15);
        EXPECT_DOUBLE_EQ(normPdf(x), normPdf(-x));
    }
}

TEST(Normal, LeftTailKeepsRelativePrecision) {
    // N(-10) = 7.619853024160527e-24. The naive (1 + erf(x/sqrt2))/2 returns 0 here.
    EXPECT_NEAR(normCdf(-10.0) / 7.619853024160527e-24, 1.0, 1e-12);
    const double naive = 0.5 * (1.0 + std::erf(-10.0 / std::sqrt(2.0)));
    EXPECT_EQ(naive, 0.0);
}

TEST(Normal, ExtremeArgumentsSaturateWithoutNaN) {
    EXPECT_EQ(normCdf(-1e300), 0.0);
    EXPECT_EQ(normCdf(1e300), 1.0);
    EXPECT_EQ(normPdf(1e300), 0.0);
    EXPECT_EQ(normCdf(-std::numeric_limits<double>::infinity()), 0.0);
    EXPECT_EQ(normCdf(std::numeric_limits<double>::infinity()), 1.0);
}
