#include "products/BarrierOption.h"
#include <algorithm>
#include <stdexcept>

namespace OptionPricer
{

    BarrierOption::BarrierOption(double K, double T, bool is_call,
                                 BarrierType barrier_type, double B,
                                 double rebate)
        : K_(K), T_(T), is_call_(is_call),
          barrier_type_(barrier_type), B_(B), rebate_(rebate)
    {
        if (K <= 0.0)
        {
            throw std::invalid_argument("Strike must be positive");
        }
        if (T <= 0.0)
        {
            throw std::invalid_argument("Maturity must be positive");
        }
        if (B <= 0.0)
        {
            throw std::invalid_argument("Barrier must be positive");
        }
        if (rebate < 0.0)
        {
            throw std::invalid_argument("Rebate must be non-negative");
        }
    }

    double BarrierOption::payoff(double S_T) const
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

    bool BarrierOption::is_breached(double S) const
    {
        switch (barrier_type_)
        {
        case BarrierType::UpAndOut:
        case BarrierType::UpAndIn:
            return S >= B_;
        case BarrierType::DownAndOut:
        case BarrierType::DownAndIn:
            return S <= B_;
        default:
            return false;
        }
    }

    const char *barrier_type_to_string(BarrierType type)
    {
        switch (type)
        {
        case BarrierType::UpAndOut:
            return "Up-and-Out";
        case BarrierType::DownAndOut:
            return "Down-and-Out";
        case BarrierType::UpAndIn:
            return "Up-and-In";
        case BarrierType::DownAndIn:
            return "Down-and-In";
        default:
            return "Unknown";
        }
    }

} // namespace OptionPricer
