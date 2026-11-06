#pragma once

namespace pricer {

/// First- and second-order sensitivities of an option price.
///
/// Conventions (all per unit change, annualised):
///  - delta = dV/dS, gamma = d2V/dS2
///  - vega  = dV/dsigma for an absolute change of 1.0 in sigma (divide by 100 for "per vol point")
///  - theta = dV/dt = -dV/dT, per year (divide by 365 for "per calendar day")
///  - rho   = dV/dr for an absolute change of 1.0 in r (divide by 100 for "per 1%")
struct Greeks {
    double delta = 0.0;
    double gamma = 0.0;
    double vega = 0.0;
    double theta = 0.0;
    double rho = 0.0;
};

}  // namespace pricer
