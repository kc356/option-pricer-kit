#include "utils/FDUtils.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace OptionPricer
{

    void solve_tridiagonal(const std::vector<double> &a,
                           const std::vector<double> &b,
                           const std::vector<double> &c,
                           const std::vector<double> &d,
                           std::vector<double> &x)
    {
        int n = static_cast<int>(b.size());
        if (n == 0)
        {
            return;
        }
        if (a.size() != static_cast<size_t>(n - 1) ||
            c.size() != static_cast<size_t>(n - 1) ||
            d.size() != static_cast<size_t>(n))
        {
            throw std::invalid_argument("Inconsistent array sizes in tridiagonal solver");
        }

        x.resize(n);

        // Thomas algorithm (LU decomposition for tridiagonal)
        std::vector<double> c_star(n - 1);
        std::vector<double> d_star(n);

        // Forward sweep
        c_star[0] = c[0] / b[0];
        d_star[0] = d[0] / b[0];

        for (int i = 1; i < n - 1; ++i)
        {
            double denom = b[i] - a[i - 1] * c_star[i - 1];
            c_star[i] = c[i] / denom;
            d_star[i] = (d[i] - a[i - 1] * d_star[i - 1]) / denom;
        }
        // Last row
        {
            int i = n - 1;
            double denom = b[i] - a[i - 1] * c_star[i - 1];
            d_star[i] = (d[i] - a[i - 1] * d_star[i - 1]) / denom;
        }

        // Back substitution
        x[n - 1] = d_star[n - 1];
        for (int i = n - 2; i >= 0; --i)
        {
            x[i] = d_star[i] - c_star[i] * x[i + 1];
        }
    }

    int solve_psor(const std::vector<double> &a,
                   const std::vector<double> &b,
                   const std::vector<double> &c,
                   const std::vector<double> &d,
                   const std::vector<double> &g,
                   std::vector<double> &x,
                   double omega,
                   double tol,
                   int max_iter)
    {
        int n = static_cast<int>(b.size());
        if (n == 0)
        {
            return 0;
        }

        // PSOR iteration
        for (int iter = 0; iter < max_iter; ++iter)
        {
            double max_diff = 0.0;

            for (int i = 0; i < n; ++i)
            {
                double x_old = x[i];

                // Gauss-Seidel update
                double sum = d[i];
                if (i > 0)
                {
                    sum -= a[i - 1] * x[i - 1];
                }
                if (i < n - 1)
                {
                    sum -= c[i] * x[i + 1];
                }

                double x_gs = sum / b[i];

                // Over-relaxation
                double x_sor = x_old + omega * (x_gs - x_old);

                // Projection (enforce x >= g for American constraint)
                x[i] = std::max(x_sor, g[i]);

                max_diff = std::max(max_diff, std::abs(x[i] - x_old));
            }

            if (max_diff < tol)
            {
                return iter + 1;
            }
        }

        return max_iter;
    }

    double linear_interp(const std::vector<double> &grid,
                         const std::vector<double> &values,
                         double x)
    {
        int n = static_cast<int>(grid.size());
        if (n == 0 || values.size() != static_cast<size_t>(n))
        {
            throw std::invalid_argument("Invalid grid or values for interpolation");
        }

        // Handle boundary cases
        if (x <= grid[0])
        {
            return values[0];
        }
        if (x >= grid[n - 1])
        {
            return values[n - 1];
        }

        // Binary search for interval
        auto it = std::lower_bound(grid.begin(), grid.end(), x);
        int j = static_cast<int>(it - grid.begin());
        if (j == 0)
        {
            j = 1;
        }

        // Linear interpolation
        double t = (x - grid[j - 1]) / (grid[j] - grid[j - 1]);
        return values[j - 1] + t * (values[j] - values[j - 1]);
    }

    std::vector<double> uniform_grid(double S_max, int N)
    {
        if (N < 2)
        {
            throw std::invalid_argument("Grid must have at least 2 points");
        }

        std::vector<double> grid(N);
        double dS = S_max / (N - 1);
        for (int i = 0; i < N; ++i)
        {
            grid[i] = i * dS;
        }
        return grid;
    }

} // namespace OptionPricer
