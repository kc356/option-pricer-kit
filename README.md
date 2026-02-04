# Option Pricer MVP

An option pricing library with C++17 computational core, Python bindings, and Streamlit UI.

## Quick Start (Docker)

Get up and running in 30 seconds:

```bash
# Build the Docker image
docker build -t option-pricer .

# Run the app
docker run -p 8501:8501 option-pricer

# Open browser: http://localhost:8501
```

## Features

### Models
- **Black-Scholes Model** with constant volatility

### Products
- **European Vanilla Options** (Call/Put)
- **Asian Arithmetic Average Options** (Call/Put)
- **American Vanilla Options** (Call/Put) with early exercise
- **Barrier Options** (Up-and-Out, Down-and-Out, Up-and-In, Down-and-In)

### Pricing Methods
- **Analytic Black-Scholes** formula for European options
- **Monte Carlo** simulation with exact lognormal discretization
  - Antithetic variates for variance reduction
  - Control variate for Asian options (using geometric Asian closed form)
- **Finite Difference (PDE)** methods
  - Crank-Nicolson scheme for European options
  - PSOR (Projected SOR) for American options with early exercise
  - Absorbing boundary conditions for Barrier options
- **Implied Volatility** calculation via bisection method

### Greeks
- Delta, Gamma, Vega, Theta, Rho for European Vanilla options (analytic)

## Architecture

The project is designed for easy extension to Heston, local volatility, and other exotic options.

```
option-pricer-kit/
├── cpp/                 # C++17 computational core
│   ├── include/         # Header files
│   │   ├── models/      # Market models
│   │   ├── products/    # Option contracts (Vanilla, Asian, American, Barrier)
│   │   ├── pricing/     # Pricing engines (Analytic, MC, FD)
│   │   ├── greeks/      # Greeks calculation
│   │   └── utils/       # Math utilities, FD solvers
│   ├── src/             # Implementation files
│   └── tests/           # Unit tests
├── python/              # pybind11 bindings
└── app/                 # Streamlit UI
```

## Build & Install

### 1. Build C++ Core

```bash
# From project root
cmake -S cpp -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# Run tests
cd build
ctest --verbose
```

### 2. Install Streamlit Dependencies

```bash
pip install -r app/requirements.txt
```

### 3. Install Python Bindings

```bash
# From project root
pip install -e python
```

This will compile the pybind11 module and install it in development mode.



## Usage

### C++ Library

```cpp
#include "models/BlackScholesModel.h"
#include "products/VanillaOption.h"
#include "pricing/AnalyticPricer.h"

BlackScholesModel model(0.05, 0.02, 0.25);  // r, q, sigma
VanillaOption call(100.0, 1.0, true);       // K, T, is_call
double price = price_bs_analytic(call, model, 100.0);
```

### Python API

```python
from option_pricer_cpp import (
    BlackScholesModel, VanillaOption, AsianOption,
    MCSettings, bs_price, bs_greeks, mc_price, 
    mc_price_asian_cv, implied_vol
)

# Price European option
model = BlackScholesModel(r=0.05, q=0.02, sigma=0.25)
option = VanillaOption(K=100.0, T=1.0, is_call=True)
price = bs_price(option, model, S0=100.0)

# Greeks
greeks = bs_greeks(option, model, S0=100.0)
print(f"Delta: {greeks['delta']:.4f}")

# Monte Carlo with antithetic variates
mc_settings = MCSettings(n_paths=100000, n_steps=252, seed=42, antithetic=True)
mc_result = mc_price(option, model, S0=100.0, settings=mc_settings)
print(f"MC Price: {mc_result['price']:.4f} ± {mc_result['std_error']:.4f}")

# Asian option with control variate
asian = AsianOption(K=100.0, T=1.0, is_call=True, n_steps=12)
asian_result = mc_price_asian_cv(asian, model, S0=100.0, settings=mc_settings, use_cv=True)

# Implied volatility
market_price = 10.45
iv = implied_vol(market_price, S0=100.0, K=100.0, T=1.0, r=0.05, q=0.02, is_call=True)
print(f"Implied Vol: {iv:.4f}")
```

### Finite Difference Pricing (Python)

```python
from option_pricer_cpp import (
    BlackScholesModel, AmericanOption, BarrierOption, BarrierType,
    FDSettings, fd_price_european, fd_price_american, fd_price_barrier
)

model = BlackScholesModel(r=0.05, q=0.02, sigma=0.25)
S0 = 100.0

# FD Settings: S_max, N (space steps), M (time steps), theta, omega, tol, max_iter
fd_settings = FDSettings(S_max=300.0, N=200, M=500, theta=0.5)

# American Put Option
american_put = AmericanOption(K=100.0, T=1.0, is_call=False)
result = fd_price_american(american_put, model, S0, fd_settings)
print(f"American Put: {result['price']:.4f}")

# Barrier Option (Down-and-Out Call)
barrier_call = BarrierOption(K=100.0, T=1.0, is_call=True, 
                              barrier=80.0, barrier_type=BarrierType.DownAndOut)
result = fd_price_barrier(barrier_call, model, S0, fd_settings)
print(f"Down-and-Out Call: {result['price']:.4f}")
```

## Finite Difference Engine

### Principle

The Finite Difference (FD) engine solves the Black-Scholes PDE numerically on a discrete grid:

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + (r-q)S\frac{\partial V}{\partial S} - rV = 0$$

We use the **Crank-Nicolson scheme** ($\theta = 0.5$), which is second-order accurate in both time and space. The scheme is unconditionally stable, making it robust for a wide range of parameter choices.

At each time step, the scheme discretizes the PDE into a tridiagonal linear system:

$$\mathbf{A} V^{n+1} = \mathbf{B} V^n + \mathbf{b}$$

where $\mathbf{A}$ and $\mathbf{B}$ are tridiagonal matrices built from the discretized differential operators, and $\mathbf{b}$ incorporates boundary conditions.

### FD Settings Guide

The `FDSettings` class controls the numerical grid and solver parameters:

| Parameter  | Description                         | Typical Values       | Guidance                                                     |
| ---------- | ----------------------------------- | -------------------- | ------------------------------------------------------------ |
| `S_max`    | Upper bound of the stock price grid | 3×K to 5×K           | Larger values improve far-OTM accuracy but reduce resolution |
| `N`        | Number of spatial steps             | 100–500              | Higher N → better accuracy, slower computation               |
| `M`        | Number of time steps                | 200–1000             | Higher M → better stability for volatile assets              |
| `theta`    | Time discretization parameter       | 0.5 (Crank-Nicolson) | 1.0 = fully implicit, 0.0 = explicit (unstable)              |
| `omega`    | PSOR relaxation parameter           | 1.0–1.5              | Only for American options; 1.0 = standard Gauss-Seidel       |
| `tol`      | PSOR convergence tolerance          | 1e-8                 | Smaller → more iterations but higher precision               |
| `max_iter` | Maximum PSOR iterations             | 1000                 | Safety limit for convergence                                 |

**Recommended starting point:**
```python
fd_settings = FDSettings(S_max=300.0, N=200, M=500, theta=0.5)
```

For higher accuracy, increase N and M proportionally. For American options, you may tune `omega` between 1.2–1.5 for faster convergence.

### American Options: PSOR Method

American options allow early exercise, introducing a **free boundary problem**. At each time step, the holder can exercise immediately for the intrinsic value or continue holding. Mathematically:

$$V(S, t) \geq \max(K - S, 0) \quad \text{(put)}$$

We use the **Projected Successive Over-Relaxation (PSOR)** algorithm to handle this constraint efficiently. At each iteration:

1. Update the option value using the linear system (standard SOR step)
2. Project the value onto the exercise constraint: $V_i = \max(V_i, \text{exercise\_value}_i)$
3. Repeat until convergence

The PSOR method converges quickly (typically 5–20 iterations per time step) and accurately captures the early exercise boundary.

### Barrier Options: Boundary Treatment

Barrier options have path-dependent features where the option becomes worthless (knock-out) or activates (knock-in) when the underlying crosses a barrier level.

**Knock-Out Options:**
- Absorbing boundary condition at the barrier: $V(\text{barrier}) = 0$
- Grid points beyond the barrier are set to zero payoff

**Knock-In Options:**
- Computed using in-out parity: $V_{\text{knock-in}} = V_{\text{vanilla}} - V_{\text{knock-out}}$
- First compute the knock-out price, then subtract from vanilla price

Supported barrier types:
- `UpAndOut` / `UpAndIn`: Barrier above current spot
- `DownAndOut` / `DownAndIn`: Barrier below current spot


### Streamlit UI

```bash
streamlit run app/streamlit_app.py
```

The UI provides:
- **Option Types**: Vanilla, Asian, American, Barrier (Up/Down-and-Out/In)
- **Pricing Methods**: Analytic (BS), Monte Carlo, Finite Difference
- **Parameters**: Interactive inputs for S0, K, T, r, q, sigma, barrier level
- **FD Settings Panel**: Configurable S_max, N, M for numerical precision
- **Greeks Display**: Delta, Gamma, Vega, Theta, Rho
- **Implied Volatility Calculator**: From market price to IV
- **Monte Carlo Convergence Plots**: Visualize price convergence
- **Payoff Diagrams**: Option value at expiry
- **FD Grid Visualization**: 3D surface of option value across (S, t)
- **Early Exercise Analysis**: American option exercise boundary plot

## Testing

The test suite validates:
- **Call-Put Parity**: European call and put prices satisfy parity relation
- **Monte Carlo Convergence**: MC prices converge to analytic prices
- **Implied Volatility**: Round-trip test (price → IV → price)
- **Control Variate Efficiency**: Asian pricing with CV reduces variance
- **FD European Pricing**: FD prices match analytic BS prices (call and put)
- **American Put Premium**: American put ≥ European put (early exercise value)
- **American Call (No Dividend)**: American call = European call when q=0
- **Barrier Knock-Out**: Barrier option price ≤ vanilla price
- **In-Out Parity**: Knock-in + Knock-out = Vanilla

```bash
cd build
ctest --verbose
```

## Numerical Stability

- Normal CDF uses `std::erfc` for stability
- Implied volatility bisection with proper bracketing [1e-6, 5.0]
- Exceptions/NaN for out-of-range inputs
- Exact lognormal path generation (not Euler discretization)
- Crank-Nicolson scheme: unconditionally stable for all parameter choices
- PSOR convergence guaranteed with proper omega (1.0–1.9)
- Thomas algorithm: O(N) tridiagonal solver with numerical stability

## Extension Points

The architecture supports easy addition of:
- **Heston Model**: Add stochastic volatility model alongside BlackScholesModel
- **Least-Squares Monte Carlo (LSM)**: Alternative method for American options
- **Double Barrier Options**: Both upper and lower barriers
- **Local Volatility**: Extend to Dupire local vol surface
- **Multi-Asset Options**: Basket, spread, and quanto options
- **Greeks for Exotics**: Extend Greeks calculation to American, Asian, Barrier options

## Project Structure Details

### C++ Core (`cpp/`)

- `include/models/` - Market models (BS, extensible to Heston)
- `include/products/` - Option contracts (Vanilla, Asian, American, Barrier)
- `include/pricing/` - Pricing engines (Analytic, MC, FD)
- `include/greeks/` - Greeks calculation
- `include/utils/` - Math utilities (CDF, implied vol, RNG, FD solvers)

### Python Bindings (`python/`)

- Uses pybind11 for zero-overhead bindings
- Automatic type conversions between C++ and NumPy
- Pythonic API with keyword arguments

### Streamlit App (`app/`)

- Clean, responsive UI
- Real-time calculations
- Interactive plots with matplotlib
- Production-ready error handling

## Performance

- C++ core optimized for speed (Release mode ~100x faster than Debug)
- Monte Carlo supports parallelization (extend with OpenMP)
- Antithetic variates reduce variance by ~50% for same computational cost
- Control variate for Asian options can reduce variance by 80-95%
- FD Thomas algorithm: O(N) per time step, total O(N×M)
- PSOR typically converges in 5-20 iterations per time step
- FD pricing: sub-second for typical grid sizes (N=200, M=500)

## License

MIT License

## References

- Hull, J. C. (2018). *Options, Futures, and Other Derivatives*
- Glasserman, P. (2003). *Monte Carlo Methods in Financial Engineering*
- [Finite Difference Methods for American Options](https://tastyhedge.com/blog/finite-difference-americans/)

## Authors

Built for quantitative finance professionals and researchers.
