#include "products/VanillaOption.h"
#include <algorithm>
#include <stdexcept>

namespace OptionPricer
{

    VanillaOption::VanillaOption(double K, double T, bool is_call)
        : K_(K), T_(T), is_call_(is_call)
    {
        if (K <= 0.0)
        {
            throw std::invalid_argument("Strike must be positive");
        }
        if (T <= 0.0)
        {
            throw std::invalid_argument("Time to maturity must be positive");
        }
    }

    double VanillaOption::payoff(double S_T) const
    {
        if (is_call_)
        {
            return std::max(S_T - K_, 0.0);
        }
        else
        {
            return std::max(K_ - S_T, 0.0);
        }
    }

} // namespace OptionPricer
