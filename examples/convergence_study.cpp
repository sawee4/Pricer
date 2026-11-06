// Monte Carlo convergence study for the Hull call (S=42, K=40, r=10%, sigma=20%, T=0.5).
//
// For each number of paths N, prices the call with M independent seeds and records:
//   - the estimate and 95% CI of the first seed (what a single run looks like),
//   - the mean reported standard error over the M runs,
//   - the root-mean-square error against the Black-Scholes price over the M runs.
// Both the SE and the RMSE should decay like 1/sqrt(N): slope -1/2 on a log-log plot.
//
// Usage: convergence_study [output.csv]   (default: mc_convergence.csv)
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "pricer/analytic_pricer.hpp"
#include "pricer/monte_carlo_pricer.hpp"

int main(int argc, char* argv[]) {
    using namespace pricer;

    const std::string path = argc > 1 ? argv[1] : "mc_convergence.csv";
    const MarketData market(42.0, 0.10, 0.0, 0.20);
    const Option call(std::make_unique<CallPayoff>(40.0), 0.5);
    const double exact = AnalyticPricer().price(call, market).price;

    const std::vector<std::size_t> pathCounts = {100,   200,   500,    1000,   2000,  5000,
                                                 10000, 20000, 50000,  100000, 200000};
    constexpr std::uint64_t kRuns = 50;

    std::ofstream csv(path);
    if (!csv) {
        std::cerr << "Cannot open " << path << '\n';
        return 1;
    }
    csv << "n_paths,price,std_error,ci_low,ci_high,bs_price,mean_std_error,rmse\n";
    csv << std::setprecision(10);

    std::cout << "Black-Scholes price: " << exact << "\n\n"
              << std::setw(8) << "N" << std::setw(12) << "price" << std::setw(12) << "SE"
              << std::setw(12) << "mean SE" << std::setw(12) << "RMSE" << '\n';

    double firstLogN = 0.0, firstLogRmse = 0.0, lastLogN = 0.0, lastLogRmse = 0.0;
    for (std::size_t i = 0; i < pathCounts.size(); ++i) {
        const std::size_t n = pathCounts[i];
        PricingResult first;
        double sumSe = 0.0, sumSqErr = 0.0;
        for (std::uint64_t seed = 1; seed <= kRuns; ++seed) {
            const PricingResult r = MonteCarloPricer(n, seed).price(call, market);
            if (seed == 1) first = r;
            sumSe += r.standardError;
            sumSqErr += (r.price - exact) * (r.price - exact);
        }
        const double meanSe = sumSe / static_cast<double>(kRuns);
        const double rmse = std::sqrt(sumSqErr / static_cast<double>(kRuns));

        csv << n << ',' << first.price << ',' << first.standardError << ',' << first.confidenceLow
            << ',' << first.confidenceHigh << ',' << exact << ',' << meanSe << ',' << rmse << '\n';
        std::cout << std::fixed << std::setprecision(5) << std::setw(8) << n << std::setw(12)
                  << first.price << std::setw(12) << first.standardError << std::setw(12)
                  << meanSe << std::setw(12) << rmse << '\n';

        const double logN = std::log(static_cast<double>(n));
        if (i == 0) {
            firstLogN = logN;
            firstLogRmse = std::log(rmse);
        }
        lastLogN = logN;
        lastLogRmse = std::log(rmse);
    }

    std::cout << "\nLog-log slope of RMSE vs N (expected -0.5): "
              << (lastLogRmse - firstLogRmse) / (lastLogN - firstLogN) << '\n'
              << "CSV written to " << path << '\n';
    return 0;
}
