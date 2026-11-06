#pragma once

#include <cstdint>
#include <random>

namespace pricer {

/// Seeded source of standard normal draws.
///
/// The engine std::mt19937_64 is fully specified by the standard, so a given seed
/// produces the same 64-bit stream everywhere. The standard *distributions*
/// (std::normal_distribution, std::uniform_real_distribution) are not: their
/// algorithms are implementation-defined, so MSVC, libstdc++ and libc++ would give
/// different prices for the same seed. We therefore map bits to (0,1) and then to
/// N(0,1) by inverse transform ourselves, which makes results bit-reproducible
/// across compilers and also keeps a monotone draw -> normal mapping
/// (convenient for antithetic variates and quasi-random sequences).
class NormalGenerator {
public:
    explicit NormalGenerator(std::uint64_t seed) : engine_(seed) {}

    /// Uniform in the open interval (0, 1), using the top 53 bits of the engine output.
    [[nodiscard]] double nextUniform();

    /// Standard normal draw.
    [[nodiscard]] double nextNormal();

private:
    std::mt19937_64 engine_;
};

}  // namespace pricer
