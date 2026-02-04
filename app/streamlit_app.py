import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import sys
import os

# Import the C++ pricing library
try:
    from option_pricer_cpp import (
        BlackScholesModel, VanillaOption, AsianOption,
        AmericanOption, BarrierOption, BarrierType,
        MCSettings, FDSettings,
        bs_price, bs_greeks, mc_price,
        mc_price_asian_cv, implied_vol,
        fd_price, fd_price_european, fd_price_american, fd_price_barrier
    )
except ImportError:
    st.error("⚠️ C++ module not found. Please build and install the package first.")
    st.code("pip install -e python")
    st.stop()

# Page configuration
st.set_page_config(
    page_title="Option Pricer MVP",
    page_icon="📈",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("📈 Option Pricer MVP")
st.markdown("*Powered by C++17 and Black-Scholes*")

# Sidebar inputs
st.sidebar.header("Option Parameters")

option_type = st.sidebar.selectbox(
    "Option Type",
    ["Vanilla", "Asian", "American", "Barrier"],
    help="European Vanilla, Asian, American, or Barrier"
)

is_call = st.sidebar.radio(
    "Call/Put",
    ["Call", "Put"],
    horizontal=True
) == "Call"

col1, col2 = st.sidebar.columns(2)
with col1:
    S0 = st.number_input("S₀ (Spot)", value=100.0, min_value=0.01, step=1.0,
                         help="Current underlying price")
    K = st.number_input("K (Strike)", value=100.0, min_value=0.01, step=1.0,
                        help="Strike price")
    T = st.number_input("T (Maturity)", value=1.0, min_value=0.01, step=0.1,
                        help="Time to maturity (years)")

with col2:
    r = st.number_input("r (Rate)", value=0.05, min_value=-0.5, max_value=0.5, step=0.01, format="%.4f",
                        help="Risk-free rate")
    q = st.number_input("q (Dividend)", value=0.02, min_value=0.0, max_value=0.5, step=0.01, format="%.4f",
                        help="Dividend yield")
    sigma = st.number_input("σ (Volatility)", value=0.25, min_value=0.001, max_value=2.0, step=0.01, format="%.4f",
                            help="Volatility (annualized)")

# Asian-specific parameters
n_monitoring = 12
if option_type == "Asian":
    n_monitoring = st.sidebar.number_input(
        "Monitoring Points",
        value=12,
        min_value=2,
        max_value=252,
        step=1,
        help="Number of monitoring points for Asian option"
    )

# Barrier-specific parameters
barrier_type_str = "Up-and-Out"
barrier_level = 120.0
barrier_type = BarrierType.UpAndOut
if option_type == "Barrier":
    barrier_type_str = st.sidebar.selectbox(
        "Barrier Type",
        ["Up-and-Out", "Down-and-Out", "Up-and-In", "Down-and-In"],
        help="Type of barrier option"
    )
    barrier_type_map = {
        "Up-and-Out": BarrierType.UpAndOut,
        "Down-and-Out": BarrierType.DownAndOut,
        "Up-and-In": BarrierType.UpAndIn,
        "Down-and-In": BarrierType.DownAndIn
    }
    barrier_type = barrier_type_map[barrier_type_str]
    
    # Default barrier level based on type
    default_barrier = 120.0 if "Up" in barrier_type_str else 80.0
    barrier_level = st.sidebar.number_input(
        "Barrier Level (B)",
        value=default_barrier,
        min_value=0.01,
        step=1.0,
        help="Barrier level that triggers knock-in/knock-out"
    )

st.sidebar.markdown("---")
st.sidebar.header("Pricing Method")

# Method selection based on option type
use_cv = True
if option_type == "Vanilla":
    method = st.sidebar.radio(
        "Method",
        ["Analytic (Black-Scholes)", "Monte Carlo", "Finite Difference"],
        help="Pricing method for vanilla options"
    )
elif option_type == "American":
    st.sidebar.info("American options priced using FD with PSOR")
    method = "Finite Difference"
elif option_type == "Barrier":
    st.sidebar.info("Barrier options priced using Finite Difference")
    method = "Finite Difference"
else:  # Asian
    method = "Monte Carlo"
    use_cv = st.sidebar.checkbox(
        "Use Control Variate",
        value=True,
        help="Use geometric Asian as control variate for variance reduction"
    )

# Finite Difference settings
fd_smax = 300.0
fd_n = 200
fd_m = 200
fd_omega = 1.2
fd_tol = 1e-6

if "Finite Difference" in method:
    st.sidebar.markdown("### Finite Difference Settings")
    fd_smax = st.sidebar.number_input(
        "S_max",
        value=300.0,
        min_value=K * 1.5,
        step=50.0,
        help="Upper bound of spot price domain (should be 3-5x strike)"
    )
    fd_n = st.sidebar.number_input(
        "N (space steps)",
        value=200,
        min_value=50,
        max_value=1000,
        step=50,
        help="Number of grid points in S direction"
    )
    fd_m = st.sidebar.number_input(
        "M (time steps)",
        value=200,
        min_value=50,
        max_value=1000,
        step=50,
        help="Number of time steps"
    )
    
    if option_type == "American":
        st.sidebar.markdown("#### PSOR Settings")
        fd_omega = st.sidebar.slider(
            "ω (relaxation)",
            min_value=1.0,
            max_value=1.9,
            value=1.2,
            step=0.05,
            help="PSOR relaxation parameter (1.0-1.8 typical)"
        )
        fd_tol = st.sidebar.number_input(
            "Tolerance",
            value=1e-6,
            format="%.1e",
            help="PSOR convergence tolerance"
        )
    
    st.sidebar.warning("⚠️ Coarse grid (small N/M) may give inaccurate results")

# Monte Carlo settings
mc_paths = 100000
mc_steps = 252
mc_seed = 42
antithetic = True

if "Monte Carlo" in method or option_type == "Asian":
    st.sidebar.markdown("### Monte Carlo Settings")
    mc_paths = st.sidebar.number_input(
        "Paths",
        value=100000,
        min_value=1000,
        max_value=10000000,
        step=10000,
        help="Number of simulation paths"
    )
    mc_steps = st.sidebar.number_input(
        "Steps",
        value=252,
        min_value=10,
        max_value=1000,
        step=10,
        help="Time steps per path"
    )
    mc_seed = st.sidebar.number_input(
        "Seed",
        value=42,
        min_value=0,
        max_value=99999,
        step=1,
        help="Random seed for reproducibility"
    )
    antithetic = st.sidebar.checkbox(
        "Antithetic Variates",
        value=True,
        help="Use antithetic variates for variance reduction"
    )

# Main area - tabs based on method
if "Finite Difference" in method:
    tab1, tab2, tab3, tab4, tab5 = st.tabs(["💰 Pricing", "📉 FD Grid", "🎯 Greeks", "📊 Implied Vol", "📈 Plots"])
else:
    tab1, tab3, tab4, tab5 = st.tabs(["💰 Pricing", "🎯 Greeks", "📊 Implied Vol", "📈 Plots"])
    tab2 = None  # No FD tab

# Create model and option
try:
    model = BlackScholesModel(r=r, q=q, sigma=sigma)
    
    if option_type == "Vanilla":
        option = VanillaOption(K=K, T=T, is_call=is_call)
    elif option_type == "Asian":
        option = AsianOption(K=K, T=T, is_call=is_call, n_steps=n_monitoring)
    elif option_type == "American":
        option = AmericanOption(K=K, T=T, is_call=is_call)
    elif option_type == "Barrier":
        option = BarrierOption(K=K, T=T, is_call=is_call, barrier_type=barrier_type, B=barrier_level)
    
    # Create FD settings if needed
    fd_settings = None
    if "Finite Difference" in method:
        fd_settings = FDSettings(
            S_max=fd_smax,
            N=fd_n,
            M=fd_m,
            theta=0.5,  # Crank-Nicolson
            omega=fd_omega,
            tol=fd_tol,
            max_iter=1000
        )
    
    # Store FD result for reuse across tabs
    fd_result = None
    
    # Tab 1: Pricing
    with tab1:
        st.header("Option Price")
        
        col_a, col_b, col_c = st.columns([2, 2, 3])
        
        with col_a:
            opt_type_display = f"{option_type} {'Call' if is_call else 'Put'}"
            if option_type == "Barrier":
                opt_type_display += f" ({barrier_type_str})"
            st.metric("Option Type", opt_type_display)
            st.metric("Spot Price", f"${S0:.2f}")
            st.metric("Strike", f"${K:.2f}")
        
        with col_b:
            st.metric("Maturity", f"{T:.2f} years")
            st.metric("Volatility", f"{sigma*100:.2f}%")
            moneyness = S0 / K
            st.metric("Moneyness (S/K)", f"{moneyness:.4f}")
            if option_type == "Barrier":
                st.metric("Barrier Level", f"${barrier_level:.2f}")
        
        # Calculate price based on option type and method
        if option_type == "Vanilla" and "Analytic" in method:
            price = bs_price(option, model, S0)
            
            with col_c:
                st.markdown("### Analytic Price")
                st.markdown(f"## ${price:.4f}")
                st.caption("Black-Scholes closed-form solution")
        
        elif option_type == "Vanilla" and "Finite Difference" in method:
            with st.spinner("Running Finite Difference..."):
                fd_result = fd_price_european(option, model, S0, fd_settings)
            
            analytic_price = bs_price(option, model, S0)
            
            with col_c:
                st.markdown("### FD Price (Crank-Nicolson)")
                st.markdown(f"## ${fd_result.price:.4f}")
                diff = abs(fd_result.price - analytic_price)
                st.caption(f"Grid: {fd_n}×{fd_m}, S_max={fd_smax}")
                st.info(f"ℹ️ Analytic price: ${analytic_price:.4f} (diff: ${diff:.6f})")
        
        elif option_type == "Vanilla" and "Monte Carlo" in method:
            mc_settings = MCSettings(
                n_paths=mc_paths,
                n_steps=mc_steps,
                seed=mc_seed,
                antithetic=antithetic
            )
            
            with st.spinner("Running Monte Carlo simulation..."):
                result = mc_price(option, model, S0, mc_settings)
            
            with col_c:
                st.markdown("### Monte Carlo Price")
                st.markdown(f"## ${result.price:.4f}")
                st.caption(f"Standard Error: ±${result.std_error:.4f}")
                st.caption(f"95% CI: [${result.ci_lower:.4f}, ${result.ci_upper:.4f}]")
                
                # Compare with analytic
                analytic_price = bs_price(option, model, S0)
                diff = abs(result.price - analytic_price)
                st.info(f"ℹ️ Analytic price: ${analytic_price:.4f} (diff: ${diff:.4f})")
        
        elif option_type == "American":
            with st.spinner("Running Finite Difference with PSOR..."):
                fd_result = fd_price_american(option, model, S0, fd_settings)
            
            # Compare with European
            euro_option = VanillaOption(K=K, T=T, is_call=is_call)
            euro_price = bs_price(euro_option, model, S0)
            early_exercise_premium = fd_result.price - euro_price
            
            with col_c:
                st.markdown("### American Price (FD + PSOR)")
                st.markdown(f"## ${fd_result.price:.4f}")
                st.caption(f"Grid: {fd_n}×{fd_m}, ω={fd_omega}")
                st.info(f"ℹ️ European equivalent: ${euro_price:.4f}")
                if early_exercise_premium > 0.001:
                    st.success(f"✅ Early exercise premium: ${early_exercise_premium:.4f}")
                else:
                    st.caption("Early exercise premium ≈ 0 (optimal to hold)")
        
        elif option_type == "Barrier":
            with st.spinner("Running Finite Difference for Barrier..."):
                fd_result = fd_price_barrier(option, model, S0, fd_settings)
            
            # Compare with vanilla
            vanilla_option = VanillaOption(K=K, T=T, is_call=is_call)
            vanilla_price = bs_price(vanilla_option, model, S0)
            
            with col_c:
                st.markdown(f"### {barrier_type_str} Price (FD)")
                st.markdown(f"## ${fd_result.price:.4f}")
                st.caption(f"Grid: {fd_n}×{fd_m}, B={barrier_level}")
                st.info(f"ℹ️ Vanilla price: ${vanilla_price:.4f}")
                
                if "Out" in barrier_type_str:
                    discount = vanilla_price - fd_result.price
                    st.caption(f"Knock-out discount: ${discount:.4f}")
                else:  # In
                    # Show in-out parity
                    out_type = BarrierType.UpAndOut if "Up" in barrier_type_str else BarrierType.DownAndOut
                    out_option = BarrierOption(K=K, T=T, is_call=is_call, barrier_type=out_type, B=barrier_level)
                    out_result = fd_price_barrier(out_option, model, S0, fd_settings)
                    st.caption(f"Out price: ${out_result.price:.4f} | Parity: In + Out ≈ Vanilla")
        
        else:  # Asian option
            mc_settings = MCSettings(
                n_paths=mc_paths,
                n_steps=n_monitoring,
                seed=mc_seed,
                antithetic=antithetic
            )
            
            with st.spinner("Running Monte Carlo simulation..."):
                if use_cv:
                    result = mc_price_asian_cv(option, model, S0, mc_settings, use_cv=True)
                    method_text = "with Control Variate"
                else:
                    result = mc_price_asian_cv(option, model, S0, mc_settings, use_cv=False)
                    method_text = "without Control Variate"
            
            with col_c:
                st.markdown("### Asian Price")
                st.markdown(f"## ${result.price:.4f}")
                st.caption(f"Standard Error: ±${result.std_error:.4f}")
                st.caption(f"95% CI: [${result.ci_lower:.4f}, ${result.ci_upper:.4f}]")
                st.caption(f"Method: {method_text}")
        
        # Additional info
        st.markdown("---")
        col_info1, col_info2, col_info3 = st.columns(3)
        
        with col_info1:
            st.markdown("**Model Parameters**")
            st.write(f"• Risk-free rate: {r*100:.2f}%")
            st.write(f"• Dividend yield: {q*100:.2f}%")
            st.write(f"• Volatility: {sigma*100:.2f}%")
        
        with col_info2:
            st.markdown("**Option Specifications**")
            st.write(f"• Type: {option_type} {'Call' if is_call else 'Put'}")
            st.write(f"• Strike: ${K:.2f}")
            st.write(f"• Maturity: {T:.2f} years")
            if option_type == "Asian":
                st.write(f"• Monitoring: {n_monitoring} points")
            if option_type == "Barrier":
                st.write(f"• Barrier: ${barrier_level:.2f} ({barrier_type_str})")
        
        with col_info3:
            st.markdown("**Pricing Details**")
            if "Monte Carlo" in method or option_type == "Asian":
                st.write(f"• Paths: {mc_paths:,}")
                st.write(f"• Steps: {mc_steps if option_type != 'Asian' else n_monitoring}")
                st.write(f"• Antithetic: {'Yes' if antithetic else 'No'}")
                if option_type == "Asian":
                    st.write(f"• Control Variate: {'Yes' if use_cv else 'No'}")
            elif "Finite Difference" in method:
                st.write(f"• Method: Crank-Nicolson")
                st.write(f"• Grid: {fd_n} × {fd_m}")
                st.write(f"• S_max: ${fd_smax:.0f}")
                if option_type == "American":
                    st.write(f"• PSOR ω: {fd_omega}")
            else:
                st.write(f"• Method: Analytic BS")
    
    # Tab 2: FD Grid (only shown for FD method)
    if tab2 is not None and "Finite Difference" in method:
        with tab2:
            st.header("Finite Difference Grid Visualization")
            
            # Get the FD result if not already computed
            if fd_result is None:
                if option_type == "Vanilla":
                    fd_result = fd_price_european(option, model, S0, fd_settings)
                elif option_type == "American":
                    fd_result = fd_price_american(option, model, S0, fd_settings)
                elif option_type == "Barrier":
                    fd_result = fd_price_barrier(option, model, S0, fd_settings)
            
            # Plot V vs S
            fig, ax = plt.subplots(figsize=(10, 6))
            
            S_grid = np.array(fd_result.grid_S)
            V_grid = np.array(fd_result.grid_V)
            payoff_grid = np.array(fd_result.grid_payoff)
            
            ax.plot(S_grid, V_grid, 'b-', linewidth=2, label='Option Value V(S,0)')
            ax.plot(S_grid, payoff_grid, 'r--', linewidth=1.5, label='Payoff at Maturity')
            
            # Mark current spot
            ax.axvline(x=S0, color='orange', linestyle='--', alpha=0.7, label=f'S₀ = ${S0:.0f}')
            ax.axvline(x=K, color='green', linestyle=':', alpha=0.7, label=f'K = ${K:.0f}')
            
            if option_type == "Barrier":
                ax.axvline(x=barrier_level, color='purple', linestyle='-', linewidth=2, alpha=0.7, 
                          label=f'Barrier B = ${barrier_level:.0f}')
            
            # Mark interpolated price
            ax.scatter([S0], [fd_result.price], color='red', s=100, zorder=5, 
                      label=f'Price at S₀: ${fd_result.price:.4f}')
            
            ax.set_xlabel('Spot Price S ($)', fontsize=12)
            ax.set_ylabel('Option Value ($)', fontsize=12)
            ax.set_title(f'{option_type} {"Call" if is_call else "Put"} - FD Solution at t=0', 
                        fontsize=14, fontweight='bold')
            ax.legend(loc='best')
            ax.grid(True, alpha=0.3)
            ax.set_xlim([0, min(fd_smax, K * 2.5)])
            ax.set_ylim([0, max(V_grid) * 1.1 if max(V_grid) > 0 else 1])
            
            st.pyplot(fig)
            plt.close()
            
            # For American options, show early exercise region
            if option_type == "American":
                st.subheader("Early Exercise Analysis")
                
                # Find early exercise region
                exercise_values = np.array([option.exercise_value(s) for s in S_grid])
                early_exercise = V_grid <= exercise_values + 0.001  # Small tolerance
                
                fig2, ax2 = plt.subplots(figsize=(10, 6))
                ax2.plot(S_grid, V_grid, 'b-', linewidth=2, label='Continuation Value')
                ax2.plot(S_grid, exercise_values, 'r--', linewidth=2, label='Immediate Exercise')
                ax2.fill_between(S_grid, 0, max(V_grid) * 1.1 if max(V_grid) > 0 else 1, 
                                where=early_exercise, alpha=0.2, color='red',
                                label='Early Exercise Region')
                
                ax2.axvline(x=S0, color='orange', linestyle='--', alpha=0.7, label=f'S₀ = ${S0:.0f}')
                ax2.set_xlabel('Spot Price S ($)', fontsize=12)
                ax2.set_ylabel('Value ($)', fontsize=12)
                ax2.set_title(f'American {"Call" if is_call else "Put"} - Early Exercise Region', 
                            fontsize=14, fontweight='bold')
                ax2.legend(loc='best')
                ax2.grid(True, alpha=0.3)
                ax2.set_xlim([0, min(fd_smax, K * 2)])
                
                st.pyplot(fig2)
                plt.close()
    
    # Tab 3: Greeks
    with tab3:
        st.header("Option Greeks")
        
        if option_type == "Vanilla":
            greeks = bs_greeks(option, model, S0)
            
            gcol1, gcol2, gcol3, gcol4, gcol5 = st.columns(5)
            
            with gcol1:
                st.metric(
                    "Delta (Δ)",
                    f"{greeks.delta:.4f}",
                    help="Sensitivity to underlying price change"
                )
                st.caption("∂V/∂S: Change in option value per $1 change in spot")
            
            with gcol2:
                st.metric(
                    "Gamma (Γ)",
                    f"{greeks.gamma:.6f}",
                    help="Curvature of Delta"
                )
                st.caption("∂²V/∂S²: Change in Delta per $1 change in spot")
            
            with gcol3:
                st.metric(
                    "Vega (ν)",
                    f"{greeks.vega:.4f}",
                    help="Sensitivity to volatility change"
                )
                st.caption("∂V/∂σ: Change in value per 1% change in volatility")
            
            with gcol4:
                st.metric(
                    "Theta (Θ)",
                    f"{greeks.theta:.4f}",
                    help="Sensitivity to time decay"
                )
                st.caption("∂V/∂t: Change in value per day of time decay")
                
            with gcol5:
                st.metric(
                    "Rho (ρ)",
                    f"{greeks.rho:.4f}",
                    help="Sensitivity to interest rate change"
                )
                st.caption("∂V/∂r: Change in value per 1% change in interest rate")
                
        else:
            st.info(f"ℹ️ Greeks for {option_type} options require numerical methods. Not implemented in current version.")
            st.markdown("**Extension Point**: Implement finite-difference Greeks (bump-and-revalue)")
    
    # Tab 4: Implied Volatility
    with tab4:
        st.header("Implied Volatility Calculator")
        
        if option_type == "Vanilla":
            ivcol1, ivcol2 = st.columns([1, 2])
            
            with ivcol1:
                st.markdown("### Input")
                market_price_input = st.number_input(
                    "Market Price",
                    value=float(bs_price(option, model, S0)),
                    min_value=0.01,
                    step=0.01,
                    format="%.4f",
                    help="Observed market price"
                )
                
                calc_iv = st.button("Calculate Implied Vol", type="primary")
            
            with ivcol2:
                st.markdown("### Result")
                
                if calc_iv:
                    try:
                        with st.spinner("Calculating implied volatility..."):
                            iv_result = implied_vol(
                                market_price_input,
                                S0, K, T, r, q, is_call
                            )
                        
                        ivcol_a, ivcol_b = st.columns(2)
                        with ivcol_a:
                            st.metric("Implied Volatility", f"{iv_result*100:.2f}%")
                            st.metric("Market Price", f"${market_price_input:.4f}")
                        
                        with ivcol_b:
                            # Verify round-trip
                            model_iv = BlackScholesModel(r=r, q=q, sigma=iv_result)
                            price_check = bs_price(option, model_iv, S0)
                            st.metric("Verification Price", f"${price_check:.4f}")
                            diff = abs(price_check - market_price_input)
                            st.metric("Difference", f"${diff:.6f}")
                        
                    except Exception as e:
                        st.error(f"❌ Error calculating implied vol: {str(e)}")
        
        else:
            st.info(f"ℹ️ Implied volatility for {option_type} options requires iterative pricing. Not implemented.")
    
    # Tab 5: Plots
    with tab5:
        st.header("Visualizations")
        
        plot_type = st.radio(
            "Select Plot",
            ["Payoff Diagram", "Monte Carlo Convergence", "Price Surface"],
            horizontal=True
        )
        
        if plot_type == "Payoff Diagram":
            st.subheader("Payoff at Maturity")
            
            # Generate spot price range
            S_range = np.linspace(K * 0.5, K * 1.5, 100)
            
            if option_type == "Vanilla":
                payoffs = np.array([option.payoff(S) for S in S_range])
                values = np.array([bs_price(option, model, S) for S in S_range])
            elif option_type == "American":
                payoffs = np.array([option.payoff(S) for S in S_range])
                values = payoffs * 1.02  # Approximate (American >= European)
            elif option_type == "Barrier":
                payoffs = np.array([option.payoff(S) for S in S_range])
                vanilla_opt = VanillaOption(K=K, T=T, is_call=is_call)
                values = np.array([bs_price(vanilla_opt, model, S) for S in S_range]) * 0.8
            else:
                # For Asian, approximate with vanilla for visualization
                st.caption("Note: Approximate payoff shown (actual depends on path)")
                payoffs = np.maximum(S_range - K, 0) if is_call else np.maximum(K - S_range, 0)
                values = payoffs * 0.8  # Rough approximation
            
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(S_range, payoffs, 'r--', label='Payoff at Maturity', linewidth=2)
            ax.plot(S_range, values, 'b-', label='Current Value (approx)', linewidth=2)
            ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
            ax.axvline(x=K, color='g', linestyle='--', alpha=0.5, label=f'Strike: ${K:.0f}')
            ax.axvline(x=S0, color='orange', linestyle='--', alpha=0.5, label=f'Current Spot: ${S0:.0f}')
            
            if option_type == "Barrier":
                ax.axvline(x=barrier_level, color='purple', linestyle='-', linewidth=2, 
                          alpha=0.7, label=f'Barrier: ${barrier_level:.0f}')
            
            ax.set_xlabel('Spot Price at Maturity ($)', fontsize=12)
            ax.set_ylabel('Value ($)', fontsize=12)
            ax.set_title(f'{option_type} {"Call" if is_call else "Put"} Payoff Diagram', fontsize=14, fontweight='bold')
            ax.legend(loc='best')
            ax.grid(True, alpha=0.3)
            
            st.pyplot(fig)
            plt.close()
        
        elif plot_type == "Monte Carlo Convergence":
            st.subheader("MC Price Convergence")
            
            if option_type not in ["Vanilla", "Asian"]:
                st.info(f"Monte Carlo not the primary method for {option_type} options. Use FD Grid tab instead.")
            else:
                if option_type == "Vanilla":
                    analytic_price = bs_price(option, model, S0)
                else:
                    analytic_price = None
                
                # Generate convergence data
                path_samples = np.logspace(3, np.log10(min(mc_paths, 500000)), 15).astype(int)
                prices = []
                errors = []
                
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                for i, n_paths in enumerate(path_samples):
                    status_text.text(f"Simulating {n_paths:,} paths...")
                    mc_settings_conv = MCSettings(
                        n_paths=n_paths,
                        n_steps=mc_steps if option_type == "Vanilla" else n_monitoring,
                        seed=mc_seed,
                        antithetic=antithetic
                    )
                    
                    if option_type == "Vanilla":
                        result = mc_price(option, model, S0, mc_settings_conv)
                    else:
                        result = mc_price_asian_cv(option, model, S0, mc_settings_conv, use_cv=use_cv)
                    
                    prices.append(result.price)
                    errors.append(result.std_error)
                    progress_bar.progress((i + 1) / len(path_samples))
                
                status_text.empty()
                progress_bar.empty()
                
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
                
                # Price convergence
                ax1.plot(path_samples, prices, 'b-', marker='o', label='MC Estimate')
                ax1.fill_between(path_samples,
                                 np.array(prices) - np.array(errors),
                                 np.array(prices) + np.array(errors),
                                 alpha=0.3, label='±1 Std Error')
                
                if analytic_price is not None:
                    ax1.axhline(y=analytic_price, color='r', linestyle='--', 
                               label=f'Analytic: ${analytic_price:.4f}')
                
                ax1.set_xlabel('Number of Paths', fontsize=12)
                ax1.set_ylabel('Estimated Price ($)', fontsize=12)
                ax1.set_title('Monte Carlo Convergence', fontsize=14, fontweight='bold')
                ax1.set_xscale('log')
                ax1.legend()
                ax1.grid(True, alpha=0.3)
                
                # Standard error decay
                ax2.plot(path_samples, errors, 'g-', marker='s')
                ax2.plot(path_samples, errors[0] * np.sqrt(path_samples[0] / path_samples), 
                        'r--', label='O(1/√n) theoretical')
                ax2.set_xlabel('Number of Paths', fontsize=12)
                ax2.set_ylabel('Standard Error ($)', fontsize=12)
                ax2.set_title('Standard Error Decay', fontsize=14, fontweight='bold')
                ax2.set_xscale('log')
                ax2.set_yscale('log')
                ax2.legend()
                ax2.grid(True, alpha=0.3)
                
                plt.tight_layout()
                st.pyplot(fig)
                plt.close()
        
        else:  # Price Surface
            st.subheader("Option Price Surface")
            
            if option_type == "Vanilla":
                # Generate grid
                S_range = np.linspace(K * 0.5, K * 1.5, 30)
                vol_range = np.linspace(0.1, 0.5, 30)
                
                S_grid_surf, vol_grid_surf = np.meshgrid(S_range, vol_range)
                price_grid = np.zeros_like(S_grid_surf)
                
                status_text = st.empty()
                status_text.text("Generating price surface...")
                
                for i in range(len(vol_range)):
                    for j in range(len(S_range)):
                        model_temp = BlackScholesModel(r=r, q=q, sigma=vol_range[i])
                        price_grid[i, j] = bs_price(option, model_temp, S_range[j])
                
                status_text.empty()
                
                fig = plt.figure(figsize=(12, 8))
                ax = fig.add_subplot(111, projection='3d')
                
                surf = ax.plot_surface(S_grid_surf, vol_grid_surf * 100, price_grid,
                                      cmap='viridis', alpha=0.8)
                
                # Mark current point
                current_price = bs_price(option, model, S0)
                ax.scatter([S0], [sigma * 100], [current_price],
                          color='r', s=100, marker='o',
                          label=f'Current: ${current_price:.2f}')
                
                ax.set_xlabel('Spot Price ($)', fontsize=11)
                ax.set_ylabel('Volatility (%)', fontsize=11)
                ax.set_zlabel('Option Price ($)', fontsize=11)
                ax.set_title(f'{option_type} {"Call" if is_call else "Put"} Price Surface',
                            fontsize=14, fontweight='bold')
                
                fig.colorbar(surf, ax=ax, shrink=0.5)
                ax.legend()
                
                st.pyplot(fig)
                plt.close()
            else:
                st.info(f"ℹ️ Price surface for {option_type} options requires extensive computation. Try the payoff or FD Grid plots instead.")

except Exception as e:
    st.error(f"❌ Error: {str(e)}")
    st.exception(e)

# Footer
st.sidebar.markdown("---")
st.sidebar.markdown("### About")
st.sidebar.info(
    """
    **Option Pricer MVP**
    
    Built with:
    - C++17 core
    - pybind11 bindings
    - Streamlit UI
    
    Features:
    - Black-Scholes pricing
    - Monte Carlo simulation
    - **Finite Difference (PDE)**
    - **American options (PSOR)**
    - **Barrier options**
    - Greeks calculation
    - Implied volatility
    - Variance reduction techniques
    """
)
