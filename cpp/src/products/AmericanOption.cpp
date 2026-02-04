#include "products/AmericanOption.h"
#include <algorithm>
#include <stdexcept>

namespace OptionPricer
{

    AmericanOption::AmericanOption(double K, double T, bool is_call)
        : K_(K), T_(T), is_call_(is_call)
    {
        if (K <= 0.0)
        {
            throw std::invalid_argument("Strike must be positive");
        }
        if (T <= 0.0)
        {
            throw std::invalid_argument("Maturity must be positive");
        }
    }

    double AmericanOption::payoff(double S) const
    {
        if (is_call_)
        {
            return std::max(S - K_, 0.0);
        }
        else
        {
            return std::max(K_ - S, 0.0);
        }
    }

    double AmericanOption::exercise_value(double S) const
    {
        // For American options, exercise value is the same as payoff
        return payoff(S);
    }

} // namespace OptionPricer
