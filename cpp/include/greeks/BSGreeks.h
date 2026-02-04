#pragma once

#include "models/BlackScholesModel.h"
#include "products/VanillaOption.h"
#include <map>
#include <string>

namespace OptionPricer
{

    /**
     * @brief Greek values for an option
     */
    struct Greeks
    {
        double delta; // ∂V/∂S
        double gamma; // ∂²V/∂S²
        double vega;  // ∂V/∂σ
        double theta; // ∂V/∂t
        double rho;   // ∂V/∂r
    };

    /**
     * @brief Calculate analytic Greeks for European vanilla options
     * @param option Vanilla option
     * @param model Black-Scholes model parameters
     * @param S0 Initial underlying price
     * @return Greeks structure
     */
    Greeks calculate_bs_greeks(const VanillaOption &option,
                               const BlackScholesModel &model,
                               double S0);

} // namespace OptionPricer
