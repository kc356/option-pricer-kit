#include "utils/MathUtils.h"
#include <cmath>

namespace OptionPricer
{
    namespace MathUtils
    {

        double norm_cdf(double x)
        {
            // Using complementary error function for numerical stability
            // CDF(x) = 0.5 * (1 + erf(x/sqrt(2)))
            // We use erfc for better accuracy: CDF(x) = 0.5 * erfc(-x/sqrt(2))
            return 0.5 * std::erfc(-x * M_SQRT1_2);
        }

        double norm_pdf(double x)
        {
            // Standard normal PDF: f(x) = (1/sqrt(2*pi)) * exp(-x²/2)
            static const double INV_SQRT_2PI = 0.3989422804014327; // 1/sqrt(2*pi)
            return INV_SQRT_2PI * std::exp(-0.5 * x * x);
        }

        double box_muller(double u1, double u2)
        {
            // Box-Muller transform
            // Z = sqrt(-2 * ln(u1)) * cos(2*pi*u2)
            return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
        }

    } // namespace MathUtils
} // namespace OptionPricer
