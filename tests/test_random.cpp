#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>

#include "pricer/normal.hpp"
#include "pricer/random.hpp"
#include "pricer/statistics.hpp"

using namespace pricer;
using pricer::math::normCdf;
using pricer::math::normInvCdf;

// ---------------------------------------------------------------------------
// Inverse normal CDF
// ---------------------------------------------------------------------------
TEST(NormInvCdf, ReferenceValues) {
    // Reference values from Python's statistics.NormalDist().inv_cdf.
    EXPECT_DOUBLE_EQ(normInvCdf(0.5), 0.0);
    EXPECT_NEAR(normInvCdf(0.975), 1.9599639845400536, 1e-15);
    EXPECT_NEAR(normInvCdf(0.3), -0.5244005127080407, 1e-15);
    EXPECT_NEAR(normInvCdf(0.02425), -1.9729610513118845, 1e-14);
    EXPECT_NEAR(normInvCdf(0.001), -3.090232306167813, 1e-14);
    EXPECT_NEAR(normInvCdf(1e-10), -6.361340902404056, 1e-13);
    EXPECT_NEAR(normInvCdf(1e-20), -9.262340089798405, 1e-13);
    EXPECT_NEAR(normInvCdf(0.999999), 4.753424308817089, 1e-11);
}

TEST(NormInvCdf, IsInverseOfCdf) {
    // Only the left half is tested: for x > 0, N(x) is rounded to a double close to 1
    // (spacing 1.1e-16), so N(x) itself has lost the information needed to recover x.
    // The right half follows from the symmetry N^{-1}(1 - p) = -N^{-1}(p).
    for (double x = -30.0; x <= 0.0; x += 0.25) {
        EXPECT_NEAR(normInvCdf(normCdf(x)), x, 1e-12 * std::max(1.0, std::abs(x))) << "x = " << x;
    }
}

TEST(NormInvCdf, EdgeCases) {
    EXPECT_EQ(normInvCdf(0.0), -std::numeric_limits<double>::infinity());
    EXPECT_EQ(normInvCdf(1.0), std::numeric_limits<double>::infinity());
    EXPECT_TRUE(std::isnan(normInvCdf(-0.1)));
    EXPECT_TRUE(std::isnan(normInvCdf(1.1)));
}

// ---------------------------------------------------------------------------
// NormalGenerator
// ---------------------------------------------------------------------------
TEST(NormalGenerator, SameSeedSameStream) {
    NormalGenerator a(123), b(123), c(124);
    bool differs = false;
    for (int i = 0; i < 1000; ++i) {
        const double x = a.nextNormal();
        EXPECT_EQ(x, b.nextNormal());
        differs = differs || (x != c.nextNormal());
    }
    EXPECT_TRUE(differs);
}

TEST(NormalGenerator, UniformsStayInOpenInterval) {
    NormalGenerator g(7);
    for (int i = 0; i < 100000; ++i) {
        const double u = g.nextUniform();
        ASSERT_GT(u, 0.0);
        ASSERT_LT(u, 1.0);
    }
}

TEST(NormalGenerator, SampleMomentsMatchStandardNormal) {
    // With n = 1e6 draws the sample mean has SE 1e-3 and the variance SE ~1.4e-3;
    // tolerances are 5 standard errors.
    NormalGenerator g(2024);
    RunningStats stats;
    double fourth = 0.0;
    const int n = 1000000;
    for (int i = 0; i < n; ++i) {
        const double z = g.nextNormal();
        stats.add(z);
        fourth += z * z * z * z;
    }
    EXPECT_NEAR(stats.mean(), 0.0, 5e-3);
    EXPECT_NEAR(stats.variance(), 1.0, 7e-3);
    EXPECT_NEAR(fourth / n, 3.0, 0.03);  // kurtosis of N(0,1)
}

TEST(NormalGenerator, IsPortableAcrossCompilers) {
    // mt19937_64 is fully specified: the 10000th output for the default seed is fixed by the
    // standard. Our transform is ours, so the first normal for seed 42 must be identical on
    // every platform (value recorded once on Linux/GCC).
    std::mt19937_64 engine;
    engine.discard(9999);
    EXPECT_EQ(engine(), 9981545732273789042ULL);
    NormalGenerator g(42);
    EXPECT_DOUBLE_EQ(g.nextNormal(), 0.69080366178458441);
}

// ---------------------------------------------------------------------------
// RunningStats (Welford)
// ---------------------------------------------------------------------------
TEST(RunningStats, KnownValues) {
    RunningStats s;
    for (double x : {2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0}) s.add(x);
    EXPECT_EQ(s.count(), 8U);
    EXPECT_DOUBLE_EQ(s.mean(), 5.0);
    EXPECT_DOUBLE_EQ(s.variance(), 32.0 / 7.0);
}

TEST(RunningStats, FewerThanTwoSamples) {
    RunningStats s;
    EXPECT_EQ(s.variance(), 0.0);
    s.add(3.0);
    EXPECT_EQ(s.variance(), 0.0);
}

TEST(RunningStats, StableWithLargeOffset) {
    // Values 1e9 + {4, 7, 13, 16}: variance 30. The naive E[X^2] - E[X]^2 formula
    // loses every digit here (1e18 squared terms vs a variance of 30).
    RunningStats s;
    double sum = 0.0, sumSq = 0.0;
    for (double x : {4.0, 7.0, 13.0, 16.0}) {
        s.add(1e9 + x);
        sum += 1e9 + x;
        sumSq += (1e9 + x) * (1e9 + x);
    }
    EXPECT_NEAR(s.variance(), 30.0, 1e-6);
    const double naive = (sumSq - sum * sum / 4.0) / 3.0;
    EXPECT_GT(std::abs(naive - 30.0), 1.0);
}
