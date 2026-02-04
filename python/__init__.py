"""
Option Pricer Python Package

High-performance option pricing library with C++ core.
"""

__version__ = "1.0.0"

try:
    from .option_pricer_cpp import (
        BlackScholesModel,
        VanillaOption,
        AsianOption,
        MCSettings,
        MCResult,
        Greeks,
        bs_price,
        bs_greeks,
        mc_price,
        mc_price_asian_cv,
        geometric_asian_price,
        implied_vol,
    )
    
    __all__ = [
        'BlackScholesModel',
        'VanillaOption',
        'AsianOption',
        'MCSettings',
        'MCResult',
        'Greeks',
        'bs_price',
        'bs_greeks',
        'mc_price',
        'mc_price_asian_cv',
        'geometric_asian_price',
        'implied_vol',
    ]
except ImportError as e:
    import warnings
    warnings.warn(f"Could not import C++ module: {e}. Please build the package first.")
