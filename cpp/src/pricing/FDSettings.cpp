#include "pricing/FDSettings.h"
#include <stdexcept>

namespace OptionPricer
{

    FDSettings::FDSettings(double S_max, int N, int M, double theta,
                           double omega, double tol, int max_iter)
        : S_max_(S_max), N_(N), M_(M), theta_(theta),
          omega_(omega), tol_(tol), max_iter_(max_iter)
    {
        if (S_max <= 0.0)
        {
            throw std::invalid_argument("S_max must be positive");
        }
        if (N < 10)
        {
            throw std::invalid_argument("N must be at least 10");
        }
        if (M < 10)
        {
            throw std::invalid_argument("M must be at least 10");
        }
        if (theta < 0.0 || theta > 1.0)
        {
            throw std::invalid_argument("theta must be in [0, 1]");
        }
        if (omega <= 0.0 || omega >= 2.0)
        {
            throw std::invalid_argument("omega must be in (0, 2)");
        }
        if (tol <= 0.0)
        {
            throw std::invalid_argument("tolerance must be positive");
        }
        if (max_iter < 1)
        {
            throw std::invalid_argument("max_iter must be at least 1");
        }
    }

} // namespace OptionPricer
