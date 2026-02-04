#include "greeks/BSGreeks.h"
#include "utils/MathUtils.h"
#include <cmath>

namespace OptionPricer
{

    Greeks calculate_bs_greeks(const VanillaOption &option,
                               const BlackScholesModel &model,
                               double S0)
    {

        double K = option.get_strike();
        double T = option.get_maturity();
        double r = model.get_r();
        double q = model.get_q();
        double sigma = model.get_sigma();

        double sigma_sqrt_T = sigma * std::sqrt(T);
        double d1 = (std::log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / sigma_sqrt_T;
        double d2 = d1 - sigma_sqrt_T;

        double pdf_d1 = MathUtils::norm_pdf(d1);
        double exp_qT = std::exp(-q * T);
        double exp_rT = std::exp(-r * T);

        double term1 = -S0 * MathUtils::norm_pdf(d1) * sigma / (2.0 * std::sqrt(T));
        double term2 = r * K * exp_rT;
        double term3 = q * S0 * exp_qT;

        Greeks greeks;

        // Delta: ∂V/∂S
        // Rho: ∂V/∂r (per 1%, so divide by 100)
        if (option.is_call())
        {
            greeks.delta = exp_qT * MathUtils::norm_cdf(d1);
            greeks.theta = term1 - term2 * MathUtils::norm_cdf(d2) + term3 * MathUtils::norm_cdf(d1);
            greeks.rho = K * T * exp_rT * MathUtils::norm_cdf(d2) / 100.0;
        }
        else
        {
            greeks.delta = -exp_qT * MathUtils::norm_cdf(-d1);
            greeks.theta = term1 + term2 * MathUtils::norm_cdf(-d2) - term3 * MathUtils::norm_cdf(-d1);
            greeks.rho = -K * T * exp_rT * MathUtils::norm_cdf(-d2) / 100.0;
        }

        // Gamma: ∂²V/∂S²
        greeks.gamma = exp_qT * pdf_d1 / (S0 * sigma_sqrt_T);

        // Vega: ∂V/∂σ (per 1%, so divide by 100)
        greeks.vega = S0 * exp_qT * pdf_d1 * std::sqrt(T) / 100.0;

        return greeks;
    }

} // namespace OptionPricer
