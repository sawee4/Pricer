// Monte Carlo Greeks study for a call (S=100, K=105, r=5%, q=2%, sigma=25%, T=1).
//
// Part 1: Greeks by every method next to the closed form.
// Part 2: noise of the bump-and-revalue delta versus the bump size h, with common random
//         numbers (CRN) and with independent draws, measured as the standard deviation of
//         the estimate over 50 seeds. Independent draws blow up like 1/h; CRN stays flat.
//         Written to a CSV for plotting.
//
// Usage: mc_greeks_study [output.csv]   (default: mc_greeks_bump.csv)
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "pricer/analytic_pricer.hpp"
#include "pricer/monte_carlo_greeks.hpp"
#include "pricer/statistics.hpp"

int main(int argc, char* argv[]) {
    using namespace pricer;
    using namespace pricer::mc_greeks;

    const std::string path = argc > 1 ? argv[1] : "mc_greeks_bump.csv";
    const MarketData market(100.0, 0.05, 0.02, 0.25);
    const Option call(std::make_unique<CallPayoff>(105.0), 1.0);
    constexpr std::size_t kPaths = 100000;
    constexpr std::uint64_t kSeed = 7;

    const Greeks exact = AnalyticPricer().greeks(call, market);
    const Greeks crn = bumpAndRevalue(MonteCarloPricer(kPaths, kSeed), call, market);
    BumpSettings independent;
    independent.commonRandomNumbers = false;
    const Greeks ind = bumpAndRevalue(MonteCarloPricer(kPaths, kSeed), call, market, independent);
    const Estimate pwDelta = pathwiseDelta(call, market, kPaths, kSeed);
    const Estimate pwVega = pathwiseVega(call, market, kPaths, kSeed);
    const Estimate lrDelta = likelihoodRatioDelta(call, market, kPaths, kSeed);

    std::cout << std::fixed << std::setprecision(5) << "N = " << kPaths << " paths\n\n"
              << std::setw(22) << "" << std::setw(11) << "delta" << std::setw(11) << "gamma"
              << std::setw(11) << "vega" << std::setw(11) << "theta" << std::setw(11) << "rho"
              << '\n';
    auto row = [](const char* name, const Greeks& g) {
        std::cout << std::setw(22) << name << std::setw(11) << g.delta << std::setw(11) << g.gamma
                  << std::setw(11) << g.vega << std::setw(11) << g.theta << std::setw(11) << g.rho
                  << '\n';
    };
    row("analytic", exact);
    row("bump, CRN", crn);
    row("bump, independent", ind);
    std::cout << "\npathwise delta         " << pwDelta.value << "  (SE " << pwDelta.standardError
              << ")\npathwise vega          " << pwVega.value << "  (SE " << pwVega.standardError
              << ")\nlikelihood ratio delta " << lrDelta.value << "  (SE " << lrDelta.standardError
              << ")\n\n";

    std::ofstream csv(path);
    if (!csv) {
        std::cerr << "Cannot open " << path << '\n';
        return 1;
    }
    csv << "relative_bump,std_delta_crn,std_delta_independent,mean_delta_crn,analytic_delta\n"
        << std::setprecision(10);
    std::cout << "Std of bump delta over 50 seeds (N = 20000 paths each)\n"
              << std::setw(10) << "h / S" << std::setw(14) << "CRN" << std::setw(14)
              << "independent" << '\n';
    for (double rel : {0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001}) {
        BumpSettings crnSettings, indSettings;
        crnSettings.relativeSpotBump = rel;
        indSettings.relativeSpotBump = rel;
        indSettings.commonRandomNumbers = false;
        RunningStats withCrn, withoutCrn;
        for (std::uint64_t seed = 1; seed <= 50; ++seed) {
            const MonteCarloPricer pricer(20000, seed);
            withCrn.add(bumpAndRevalue(pricer, call, market, crnSettings).delta);
            withoutCrn.add(bumpAndRevalue(pricer, call, market, indSettings).delta);
        }
        const double stdCrn = std::sqrt(withCrn.variance());
        const double stdInd = std::sqrt(withoutCrn.variance());
        csv << rel << ',' << stdCrn << ',' << stdInd << ',' << withCrn.mean() << ','
            << exact.delta << '\n';
        std::cout << std::setw(10) << std::setprecision(4) << rel << std::setw(14)
                  << std::setprecision(5) << stdCrn << std::setw(14) << stdInd << '\n';
    }
    std::cout << "\nCSV written to " << path << '\n';
    return 0;
}
