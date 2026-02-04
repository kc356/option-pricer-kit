#pragma once

#include <vector>

namespace OptionPricer
{

    /**
     * @brief Solve tridiagonal linear system using Thomas algorithm
     *
     * Solves the system: A * x = d
     * where A is tridiagonal with sub-diagonal a, diagonal b, super-diagonal c
     *
     * @param a Sub-diagonal (size n-1)
     * @param b Main diagonal (size n)
     * @param c Super-diagonal (size n-1)
     * @param d Right-hand side (size n)
     * @param x Solution (size n) - output
     */
    void solve_tridiagonal(const std::vector<double> &a,
                           const std::vector<double> &b,
                           const std::vector<double> &c,
                           const std::vector<double> &d,
                           std::vector<double> &x);

    /**
     * @brief PSOR (Projected Successive Over-Relaxation) solver for American options
     *
     * Solves: A * x = d subject to x >= g (exercise constraint)
     *
     * @param a Sub-diagonal (size n-1)
     * @param b Main diagonal (size n)
     * @param c Super-diagonal (size n-1)
     * @param d Right-hand side (size n)
     * @param g Lower bound (exercise value) (size n)
     * @param x Solution (size n) - input/output (initial guess on input)
     * @param omega Relaxation parameter (1.0-1.8)
     * @param tol Convergence tolerance
     * @param max_iter Maximum iterations
     * @return Number of iterations used
     */
    int solve_psor(const std::vector<double> &a,
                   const std::vector<double> &b,
                   const std::vector<double> &c,
                   const std::vector<double> &d,
                   const std::vector<double> &g,
                   std::vector<double> &x,
                   double omega,
                   double tol,
                   int max_iter);

    /**
     * @brief Linear interpolation in a 1D array
     *
     * @param grid Grid points (strictly increasing)
     * @param values Values at grid points
     * @param x Point to interpolate at
     * @return Interpolated value
     */
    double linear_interp(const std::vector<double> &grid,
                         const std::vector<double> &values,
                         double x);

    /**
     * @brief Generate uniform grid from 0 to S_max
     *
     * @param S_max Upper bound
     * @param N Number of grid points
     * @return Grid vector
     */
    std::vector<double> uniform_grid(double S_max, int N);

} // namespace OptionPricer
