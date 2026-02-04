#pragma once

#include "models/BlackScholesModel.h"
#include "products/VanillaOption.h"
#include "products/AmericanOption.h"
#include "products/BarrierOption.h"
#include "pricing/FDSettings.h"
#include <vector>

namespace OptionPricer
{

    /**
     * @brief Result from FD pricing including grid data for visualization
     */
    struct FDResult
    {
        double price;                    // Option price at S0
        std::vector<double> grid_S;      // Spot price grid
        std::vector<double> grid_V;      // Option values at t=0
        std::vector<double> grid_payoff; // Payoff values (for American comparison)
    };

    /**
     * @brief Finite Difference pricing for European vanilla option
     *
     * Uses Crank-Nicolson scheme with Thomas algorithm.
     *
     * @param option European vanilla option
     * @param model Black-Scholes model parameters
     * @param S0 Initial spot price
     * @param settings FD settings (grid, scheme parameters)
     * @return FDResult with price and grid data
     */
    FDResult fd_price_european(const VanillaOption &option,
                               const BlackScholesModel &model,
                               double S0,
                               const FDSettings &settings);

    /**
     * @brief Finite Difference pricing for American vanilla option
     *
     * Uses Crank-Nicolson scheme with PSOR for early exercise constraint.
     *
     * @param option American vanilla option
     * @param model Black-Scholes model parameters
     * @param S0 Initial spot price
     * @param settings FD settings (grid, scheme, PSOR parameters)
     * @return FDResult with price and grid data
     */
    FDResult fd_price_american(const AmericanOption &option,
                               const BlackScholesModel &model,
                               double S0,
                               const FDSettings &settings);

    /**
     * @brief Finite Difference pricing for European barrier option
     *
     * Uses Crank-Nicolson with absorbing boundary at barrier level.
     * For knock-in options, uses in-out parity: In = Vanilla - Out
     *
     * @param option Barrier option
     * @param model Black-Scholes model parameters
     * @param S0 Initial spot price
     * @param settings FD settings
     * @return FDResult with price and grid data
     */
    FDResult fd_price_barrier(const BarrierOption &option,
                              const BlackScholesModel &model,
                              double S0,
                              const FDSettings &settings);

    /**
     * @brief Convenience function for simple price extraction
     */
    double fd_price(const VanillaOption &option,
                    const BlackScholesModel &model,
                    double S0,
                    const FDSettings &settings);

    double fd_price(const AmericanOption &option,
                    const BlackScholesModel &model,
                    double S0,
                    const FDSettings &settings);

    double fd_price(const BarrierOption &option,
                    const BlackScholesModel &model,
                    double S0,
                    const FDSettings &settings);

} // namespace OptionPricer
