#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

namespace OptionPricer
{
    namespace MathUtils
    {

        /**
         * @brief Standard normal cumulative distribution function
         * @param x Input value
         * @return Probability P(Z <= x) where Z ~ N(0,1)
         */
        double norm_cdf(double x);

        /**
         * @brief Standard normal probability density function
         * @param x Input value
         * @return Density f(x) where x ~ N(0,1)
         */
        double norm_pdf(double x);

        /**
         * @brief Generate standard normal random variable (Box-Muller)
         * @param u1 Uniform random variable [0,1]
         * @param u2 Uniform random variable [0,1]
         * @return Standard normal random variable
         */
        double box_muller(double u1, double u2);

    } // namespace MathUtils
} // namespace OptionPricer
