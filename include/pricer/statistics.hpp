#pragma once

#include <cstddef>

namespace pricer {

/// One-pass mean and variance accumulator (Welford's algorithm).
///
/// The textbook formula Var = E[X^2] - E[X]^2 subtracts two large, nearly equal
/// numbers and loses all significant digits when the mean is large compared with
/// the standard deviation. Welford's update is numerically stable and needs O(1) memory.
class RunningStats {
public:
    void add(double x) noexcept {
        ++count_;
        const double delta = x - mean_;
        mean_ += delta / static_cast<double>(count_);
        m2_ += delta * (x - mean_);
    }

    [[nodiscard]] std::size_t count() const noexcept { return count_; }
    [[nodiscard]] double mean() const noexcept { return mean_; }

    /// Unbiased sample variance (divides by n - 1). Zero when fewer than two samples.
    [[nodiscard]] double variance() const noexcept {
        return count_ > 1 ? m2_ / static_cast<double>(count_ - 1) : 0.0;
    }

private:
    std::size_t count_ = 0;
    double mean_ = 0.0;
    double m2_ = 0.0;  ///< sum of squared deviations from the current mean
};

/// One-pass accumulator for the covariance of a pair (X, Y), also by Welford-style updates.
/// Used by the control-variate estimator, which needs Var(Y), Var(X) and Cov(X, Y).
class RunningCovariance {
public:
    void add(double x, double y) noexcept {
        ++count_;
        const double n = static_cast<double>(count_);
        const double dx = x - meanX_;
        meanX_ += dx / n;
        const double dy = y - meanY_;
        meanY_ += dy / n;
        m2X_ += dx * (x - meanX_);
        m2Y_ += dy * (y - meanY_);
        cXY_ += dx * (y - meanY_);  // uses the old mean of x and the new mean of y
    }

    [[nodiscard]] std::size_t count() const noexcept { return count_; }
    [[nodiscard]] double meanX() const noexcept { return meanX_; }
    [[nodiscard]] double meanY() const noexcept { return meanY_; }
    [[nodiscard]] double varianceX() const noexcept { return unbiased(m2X_); }
    [[nodiscard]] double varianceY() const noexcept { return unbiased(m2Y_); }
    [[nodiscard]] double covariance() const noexcept { return unbiased(cXY_); }

private:
    [[nodiscard]] double unbiased(double m) const noexcept {
        return count_ > 1 ? m / static_cast<double>(count_ - 1) : 0.0;
    }

    std::size_t count_ = 0;
    double meanX_ = 0.0;
    double meanY_ = 0.0;
    double m2X_ = 0.0;
    double m2Y_ = 0.0;
    double cXY_ = 0.0;
};

}  // namespace pricer
