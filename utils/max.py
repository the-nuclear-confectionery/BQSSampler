import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar


def read_momentum_grid(filepath):
    df = pd.read_csv(filepath, sep='\s+', header=None)
    return df.iloc[:, 0].values


def weight_function(p, mbar, chem):
    E = np.sqrt(p**2 + mbar**2)
    return np.exp(p) / (np.exp(E - chem) - 1)


def max_weight_condition(mbar, chem, p_grid):
    """
    Computes the maximum of w_n(p) for a given (mbar, chem), 
    but excludes cases where chem >= E (divergence).
    """
    E_values = np.sqrt(p_grid**2 + mbar**2)
    if np.any(chem >= E_values):  # Check if chem >= E for any p
        return np.nan  # Assign NaN to avoid invalid calculations
    
    w_values = np.array([weight_function(p, mbar, chem) for p in p_grid])
    return np.nanmax(w_values)  # Ensure NaN propagates correctly


def compute_max_w_grid(mbar_vals, chem_vals, p_grid):
    max_w_values = np.zeros((len(chem_vals), len(mbar_vals)))

    for i, chem in enumerate(chem_vals):
        for j, mbar in enumerate(mbar_vals):
            max_w_values[i, j] = max_weight_condition(mbar, chem, p_grid)

    return max_w_values


def fit_boundary(mbar_values, chem_values):
    mbar_values = np.array(mbar_values)
    chem_values = np.array(chem_values)

    mask_neg = chem_values < 0
    mask_pos = chem_values > 0

    poly_neg = np.polyfit(mbar_values[mask_neg], chem_values[mask_neg], deg=4)
    poly_pos = np.polyfit(mbar_values[mask_pos], chem_values[mask_pos], deg=4)

    return poly_neg, poly_pos, np.poly1d(poly_neg), np.poly1d(poly_pos)


def plot_max_w_grid(mbar_vals, chem_vals, max_w_values, mbar_chem_pairs):
    """
    Plots the colormap of max w_n(p) over the full (m/T, chem) grid,
    excluding invalid points where chem >= E.
    """
    mbar_values, chem_values = zip(*mbar_chem_pairs)
    poly_neg, poly_pos, poly_func_neg, poly_func_pos = fit_boundary(mbar_values, chem_values)

    print("Polynomial coefficients for chem < 0:", poly_neg)
    print("Polynomial coefficients for chem > 0:", poly_pos)

    plt.figure(figsize=(8, 6))
    
    # Mask invalid (NaN) values to ensure proper plotting
    max_w_values_masked = np.ma.masked_invalid(max_w_values)

    # Plot colormap
    # Define the limits for color saturation
    vmin = 0      # Minimum value for colormap
    vmax = 5      # Maximum value to avoid dominance of large values
    
    plt.figure(figsize=(8, 6))
    

    vmin = 0  # Minimum value for colormap
    vmax = 3  # Fixed maximum value, anything above this will be clamped

    # Clamp the values at vmax to ensure anything above it appears the same color
    clamped_w_values = np.clip(max_w_values, vmin, vmax)

    plt.figure(figsize=(8, 6))

    # Plot colormap with fixed limits
    contour = plt.contourf(mbar_vals, chem_vals, clamped_w_values, levels=100, cmap='viridis', vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(contour)
    cbar.set_label("Max w_n(p)")
    plt.contour(mbar_vals, chem_vals, max_w_values, levels=[1], colors='red', linewidths=2, linestyles='dashed', label="w_n(p) = 1")
    # Ensure correct color scaling
    cbar.mappable.set_clim(vmin, vmax)
    

    # Plot polynomial boundaries in red
    mbar_fit = np.linspace(min(mbar_vals), max(mbar_vals), 500)
    #chem_fit_neg = poly_func_neg(mbar_fit)
    #chem_fit_pos = poly_func_pos(mbar_fit)

    #intersection_idx = np.argmin(np.abs(chem_fit_neg - chem_fit_pos))
    #mbar_cutoff = mbar_fit[intersection_idx]

    # Valid regions only
   # valid_mbar_neg = mbar_fit[mbar_fit <= mbar_cutoff]
   # valid_mbar_pos = mbar_fit[mbar_fit >= mbar_cutoff]
    #valid_chem_neg = chem_fit_neg[:len(valid_mbar_neg)]
    #valid_chem_pos = chem_fit_pos[len(mbar_fit) - len(valid_mbar_pos):]

    #plt.plot(valid_mbar_neg, valid_chem_neg, 'r-', linewidth=2, label="Poly fit (chem < 0)")
    #plt.plot(valid_mbar_pos, valid_chem_pos, 'r-', linewidth=2, label="Poly fit (chem > 0)")

    plt.xlabel("m/T")
    plt.ylabel("chem = (mu_Q + mu_B + mu_S) / T")
    plt.title("Maximum w_n(p) over full (m/T, chem) grid")
    plt.legend()
    plt.grid(True)
    plt.savefig("max_w_grid.png")
    plt.show()


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_gauss_laguerre_grid>")
        sys.exit(1)

    filepath = sys.argv[1]
    p_grid = read_momentum_grid(filepath)

    # Define the (m/T, chem) grid
    mbar_vals = np.linspace(0.1, 7.5, 1000)
    chem_vals = np.linspace(-7, 7, 1000)

    # Compute the maximum weight function over the grid
    max_w_values = compute_max_w_grid(mbar_vals, chem_vals, p_grid)

    # Find the boundary curve
    mbar_chem_pairs = [(m, c) for i, c in enumerate(chem_vals) for j, m in enumerate(mbar_vals) if max_w_values[i, j] > 1]

    if mbar_chem_pairs:
        plot_max_w_grid(mbar_vals, chem_vals, max_w_values, mbar_chem_pairs)
    else:
        print("No valid (mbar, chem) pairs found.")
