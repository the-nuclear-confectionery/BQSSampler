import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar


def read_momentum_grid(filepath):
    df = pd.read_csv(filepath, sep='\s+', header=None)
    return df.iloc[:, 0].values

def weight_function2(p, mbar, chem):
    E = np.sqrt(p**2 + mbar**2)
    return np.exp(p) / (np.exp(E - chem) - 1)


def max_weight_condition2(mbar, chem, p_grid):
    w_values = np.array([weight_function2(p, mbar, chem) for p in p_grid])
    return np.max(w_values) - 1


def find_all_mbar_chem(filepath, tol=1e-6):
    p_grid = read_momentum_grid(filepath)
    max_p = np.max(p_grid)
    
    mbar_chem_pairs = []
    chem_values = np.linspace(-10, 10, 10000)
    if 0 not in chem_values:
        chem_values = np.append(chem_values, 0)
        chem_values = np.sort(chem_values)
    
    for chem in chem_values:
        lower_bound = max(0.1, chem - max_p)
        upper_bound = lower_bound + 9.0
        
        f_a = max_weight_condition2(lower_bound, chem, p_grid)
        f_b = max_weight_condition2(upper_bound, chem, p_grid)
        
        if f_a * f_b > 0:
            for new_b in np.linspace(upper_bound + 1.0, upper_bound + 5.0, 5):
                f_new_b = max_weight_condition2(new_b, chem, p_grid)
                if f_a * f_new_b < 0:
                    upper_bound = new_b
                    break
            else:
                continue
        
        sol = root_scalar(lambda mbar: max_weight_condition2(mbar, chem, p_grid),
                          bracket=[lower_bound, upper_bound], method='brentq', xtol=tol)
        
        if sol.converged:
            mbar_value = sol.root
            mbar_chem_pairs.append((mbar_value, chem))
    
    return mbar_chem_pairs


def fit_boundary(mbar_values, chem_values):
    mbar_values = np.array(mbar_values)
    chem_values = np.array(chem_values)
    
    mask_neg = chem_values < 0
    mask_pos = chem_values > 0
    
    poly_neg = np.polyfit(mbar_values[mask_neg], chem_values[mask_neg], deg=4)
    poly_pos = np.polyfit(mbar_values[mask_pos], chem_values[mask_pos], deg=4)
    
    return poly_neg, poly_pos, np.poly1d(poly_neg), np.poly1d(poly_pos)


def plot_mbar_chem_pairs(mbar_chem_pairs):
    mbar_values, chem_values = zip(*mbar_chem_pairs)
    poly_neg, poly_pos, poly_func_neg, poly_func_pos = fit_boundary(mbar_values, chem_values)
    
    print("Polynomial coefficients for chem < 0:", poly_neg)
    print("Polynomial coefficients for chem > 0:", poly_pos)
    plt.figure(figsize=(8, 6))
    plt.scatter(mbar_values, chem_values, s=15, alpha=0.7, label="Boundary Points")
    
    mbar_fit = np.linspace(min(mbar_values), max(mbar_values), 500)
    chem_fit_neg = poly_func_neg(mbar_fit)
    chem_fit_pos = poly_func_pos(mbar_fit)
    
    intersection_idx = np.argmin(np.abs(chem_fit_neg - chem_fit_pos))
    mbar_cutoff = mbar_fit[intersection_idx]
    
    mbar_fit_neg = mbar_fit[mbar_fit <= mbar_cutoff]
    mbar_fit_pos = mbar_fit[mbar_fit >= mbar_cutoff]
    chem_fit_neg = chem_fit_neg[:len(mbar_fit_neg)]
    chem_fit_pos = chem_fit_pos[len(mbar_fit) - len(mbar_fit_pos):]
    
    plt.plot(mbar_fit_neg, chem_fit_neg, 'r-', label="Fit (chem < 0)")
    plt.plot(mbar_fit_pos, chem_fit_pos, 'g-', label="Fit (chem > 0)")
    
    plt.fill_between(mbar_fit_neg, np.max(chem_values), chem_fit_neg, color='red', alpha=0.3, label='w_n(p) > 1')
    plt.fill_between(mbar_fit_pos, chem_fit_pos, np.max(chem_values), color='red', alpha=0.3)
    
    plt.xlabel("m/T")
    plt.ylabel("chem = (mu_Q + mu_B + mu_S) / T")
    plt.title("(m/T, chem) Pairs where max w_n(p) = 1")
    plt.legend()
    plt.grid(True)
    plt.savefig("mbar_chem_pairs2.png")
    plt.show()


def weight_function(p, mbar, chem):
    E = np.sqrt(p**2 + mbar**2)
    return np.exp(p) / (np.exp(E - chem) - 1)


def max_weight_condition(mbar, chem, p_grid):
    """
    Computes the maximum of w_n(p) for a given (mbar, chem), 
    but excludes cases where chem >= E (divergence).
    """
    E_values = np.sqrt(p_grid**2 + mbar**2)
    #if np.any(chem >= E_values):  # Check if chem >= E for any p
    #    return np.nan  # Assign NaN to avoid invalid calculations
    
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


def plot_max_w_grid(mbar_vals, chem_vals, max_w_values, mbar_chem_pairs,mbar_chem_pairs2):
    """
    Plots the colormap of max w_n(p) over the full (m/T, chem) grid,
    excluding invalid points where chem >= E.
    """
    mbar_values, chem_values = zip(*mbar_chem_pairs)

    plt.figure(figsize=(8, 6))
    
    # Mask invalid (NaN) values to ensure proper plotting
    max_w_values_masked = np.ma.masked_invalid(max_w_values)

    # Plot colormap
    # Define the limits for color saturation
    vmin = 0      # Minimum value for colormap
    vmax = 3      # Maximum value to avoid dominance of large values

    # Clamp the values at vmax to ensure anything above it appears the same color
    clamped_w_values = np.clip(max_w_values, vmin, vmax)


    # Plot colormap with fixed limits
    contour = plt.contourf(mbar_vals, chem_vals, clamped_w_values, levels=100, cmap='jet', vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(contour)
    cbar.set_label(r"Max $w_n(\Omega,m/T, p/T)$")
    plt.contour(mbar_vals, chem_vals, max_w_values, levels=[1], colors='black', linewidths=2, linestyles='dashed', label="w_n(p) = 1",zorder=10)
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
    mbar_values2, chem_values2 = zip(*mbar_chem_pairs2)
    poly_neg, poly_pos, poly_func_neg, poly_func_pos = fit_boundary(mbar_values2, chem_values2)
    
    print("Polynomial coefficients for chem < 0:", poly_neg)
    print("Polynomial coefficients for chem > 0:", poly_pos)

    #plt.scatter(mbar_values2, chem_values2, s=15, alpha=0.7, label="Boundary Points")
    
    mbar_fit = np.linspace(min(mbar_values), max(mbar_values), 500)
    chem_fit_neg = poly_func_neg(mbar_fit)
    chem_fit_pos = poly_func_pos(mbar_fit)
    
    intersection_idx = np.argmin(np.abs(chem_fit_neg - chem_fit_pos))
    mbar_cutoff = 0.8551
    
    mbar_fit_neg = mbar_fit[mbar_fit <= mbar_cutoff]
    mbar_fit_pos = mbar_fit[mbar_fit >= mbar_cutoff]
    chem_fit_neg = chem_fit_neg[:len(mbar_fit_neg)]
    chem_fit_pos = chem_fit_pos[len(mbar_fit) - len(mbar_fit_pos):]
    
    plt.plot(mbar_fit_neg, chem_fit_neg, 'darkorange', lw=3, label="Fit")
    plt.plot(mbar_fit_pos, chem_fit_pos, 'darkorange', lw=3)

    # Compute the E = chem boundary (dotted line)
    E_grid = np.sqrt(mbar_vals[:, None]**2 + p_grid[None, :]**2)
    E_min = np.min(E_grid, axis=1)  # Find min E for each mbar
    plt.fill_between(mbar_vals, E_min, np.max(chem_vals), color='white', alpha=0.4, hatch='//')
    # Plot the E = chem boundary
    plt.plot(mbar_vals, E_min, 'k:',color="white" ,lw=2, label=r"$E/T = \Omega$")  # Dotted line
    plt.scatter(0.92169, 0.0380237) 
    #plt.fill_between(mbar_fit_neg, np.max(chem_values2), chem_fit_neg, color='red', alpha=0.3, label='w_n(p) > 1')
    #plt.fill_between(mbar_fit_pos, chem_fit_pos, np.max(chem_values2), color='red', alpha=0.3)
    #plt.scatter(0.92006,0.0341833, s=10, color='red',label='Boundary Points')
    plt.xlabel(r"$m/T$")
    plt.ylabel(r"$ \Omega = (\mu_Q N_Q + \mu_B N_B + \mu_S N_S) / T$")
    plt.legend()
    plt.ylim(-2,2)
    plt.savefig("max_w_grid2.png")
    plt.show()


def save_max_w_table(mbar_vals, chem_vals, max_w_values, filename="max_w_table.dat"):
    """
    Saves a structured table of (m/T, chem, max w_n(p)) where max w_n(p) > 1,
    ensuring correct order without using flattening.
    Adds a header line with # Nx Ny Npoints.
    """
    Nx = len(mbar_vals)
    Ny = len(chem_vals)
    Npoints = Nx * Ny  # Total number of grid points

    with open(filename, 'w') as f:
        # Write the header
        f.write(f"# {Nx} {Ny} {Npoints}\n")

        # Manually iterate over indices to maintain correct structure
        for j, mbar in enumerate(mbar_vals):  # Iterate over mbar first (outer loop)
            for i, chem in enumerate(chem_vals):  # Iterate over chem second (inner loop)
                max_w = max_w_values[i, j]  # Access max_w_values correctly
                f.write(f"{mbar:.6f} {chem:.6f} {max_w:.6f}\n")  # Write line

    print(f"Filtered table saved to {filename}")



if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_gauss_laguerre_grid>")
        sys.exit(1)

    filepath = sys.argv[1]
    p_grid = read_momentum_grid(filepath)

    # Define the (m/T, chem) grid
    mbar_vals = np.linspace(0.1, 20., 1001)
    chem_vals = np.linspace(-7, 7, 1001)

    # Compute the maximum weight function over the grid
    max_w_values = compute_max_w_grid(mbar_vals, chem_vals, p_grid)

    # Find the boundary curve
    mbar_chem_pairs = [(m, c) for i, c in enumerate(chem_vals) for j, m in enumerate(mbar_vals) if max_w_values[i, j] > 1]
    mbar_chem_pairs2 = find_all_mbar_chem(filepath)

    if mbar_chem_pairs:
        plot_max_w_grid(mbar_vals, chem_vals, max_w_values, mbar_chem_pairs,mbar_chem_pairs2)
    else:
        print("No valid (mbar, chem) pairs found.")


    
    if mbar_chem_pairs2:
        print("All (mbar, chem) pairs where max w_n(p) = 1:")
        for mbar, chem in mbar_chem_pairs2:
            print(f"(mbar, chem) = ({mbar:.4f}, {chem:.4f})")
        
        #plot_mbar_chem_pairs(mbar_chem_pairs2)
    else:
        print("No valid (mbar, chem) pairs found.")

    # Save the filtered table where max_w > 1
    save_max_w_table(mbar_vals, chem_vals, max_w_values)

