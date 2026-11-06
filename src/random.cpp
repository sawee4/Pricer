#include "pricer/random.hpp"

#include "pricer/normal.hpp"

namespace pricer {

double NormalGenerator::nextUniform() {
    // (k + 0.5) / 2^53 with k in [0, 2^53): never exactly 0 or 1, so the inverse CDF is finite.
    constexpr double kInv2Pow53 = 1.0 / 9007199254740992.0;
    const std::uint64_t k = engine_() >> 11;
    return (static_cast<double>(k) + 0.5) * kInv2Pow53;
}

double NormalGenerator::nextNormal() { return math::normInvCdf(nextUniform()); }

}  // namespace pricer
