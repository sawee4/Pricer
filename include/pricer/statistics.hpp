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

}  // namespace pricer
