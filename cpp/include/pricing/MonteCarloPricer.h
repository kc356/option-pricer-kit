#pragma once

#include "models/BlackScholesModel.h"
#include "products/VanillaOption.h"
#include "products/AsianOption.h"
#include "pricing/MCSettings.h"
#include <map>
#include <string>

namespace OptionPricer
{

    /**
     * @brief Result structure for Monte Carlo pricing
     */
    struct MCResult
    {
        double price;     // Estimated price
        double std_error; // Standard error of the mean
        double ci_lower;  // 95% confidence interval lower bound
        double ci_upper;  // 95% confidence interval upper bound
    };

    /**
     * @brief Monte Carlo pricing for vanilla options
     * @param option Vanilla option to price
     * @param model Black-Scholes model parameters
     * @param S0 Initial underlying price
     * @param settings MC simulation settings
     * @return MC result with price and statistics
     */
    MCResult mc_price(const VanillaOption &option,
                      const BlackScholesModel &model,
                      double S0,
                      const MCSettings &settings);

    /**
     * @brief Monte Carlo pricing for Asian options
     * @param option Asian option to price
     * @param model Black-Scholes model parameters
     * @param S0 Initial underlying price
     * @param settings MC simulation settings
     * @return MC result with price and statistics
     */
    MCResult mc_price(const AsianOption &option,
                      const BlackScholesModel &model,
                      double S0,
                      const MCSettings &settings);

    /**
     * @brief Monte Carlo pricing for Asian options with control variate
     *
     * Uses geometric Asian option as control variate for variance reduction
     *
     * @param option Asian option to price
     * @param model Black-Scholes model parameters
     * @param S0 Initial underlying price
     * @param settings MC simulation settings
     * @param use_cv Whether to use control variate
     * @return MC result with price and statistics
     */
    MCResult mc_price_asian_cv(const AsianOption &option,
                               const BlackScholesModel &model,
                               double S0,
                               const MCSettings &settings,
                               bool use_cv = true);

} // namespace OptionPricer
