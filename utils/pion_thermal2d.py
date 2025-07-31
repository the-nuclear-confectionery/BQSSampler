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
    w_values = np.array([weight_function(p, mbar, chem) for p in p_grid])
    return np.max(w_values) - 1


def find_all_mbar_chem(filepath, tol=1e-6):
    p_grid = read_momentum_grid(filepath)
    max_p = np.max(p_grid)
    
    mbar_chem_pairs = []
    chem_values = np.linspace(-6, 6, 10000)
    if 0 not in chem_values:
        chem_values = np.append(chem_values, 0)
        chem_values = np.sort(chem_values)
    
    for chem in chem_values:
        lower_bound = max(0.1, chem - max_p)
        upper_bound = lower_bound + 2.0
        
        f_a = max_weight_condition(lower_bound, chem, p_grid)
        f_b = max_weight_condition(upper_bound, chem, p_grid)
        
        if f_a * f_b > 0:
            for new_b in np.linspace(upper_bound + 1.0, upper_bound + 5.0, 5):
                f_new_b = max_weight_condition(new_b, chem, p_grid)
                if f_a * f_new_b < 0:
                    upper_bound = new_b
                    break
            else:
                continue
        
        sol = root_scalar(lambda mbar: max_weight_condition(mbar, chem, p_grid),
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
    mbar_cutoff = 0.8551
    
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
    plt.savefig("mbar_chem_pairs.png")
    plt.show()


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_gauss_laguerre_grid>")
        sys.exit(1)
    
    filepath = sys.argv[1]
    mbar_chem_pairs = find_all_mbar_chem(filepath)
    
    if mbar_chem_pairs:
        print("All (mbar, chem) pairs where max w_n(p) = 1:")
        for mbar, chem in mbar_chem_pairs:
            print(f"(mbar, chem) = ({mbar:.4f}, {chem:.4f})")
        
        plot_mbar_chem_pairs(mbar_chem_pairs)
    else:
        print("No valid (mbar, chem) pairs found.")
