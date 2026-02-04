#include "models/BlackScholesModel.h"
#include <stdexcept>

namespace OptionPricer
{

    BlackScholesModel::BlackScholesModel(double r, double q, double sigma)
        : r_(r), q_(q), sigma_(sigma)
    {
        if (sigma <= 0.0)
        {
            throw std::invalid_argument("Volatility must be positive");
        }
    }

} // namespace OptionPricer
