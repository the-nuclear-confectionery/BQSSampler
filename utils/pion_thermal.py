import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

def read_momentum_grid(filepath):
    """ Reads the first column from the Gauss-Laguerre grid file """
    df = pd.read_csv(filepath, sep='\s+', header=None)
    return df.iloc[:, 0].values  # Extract only the first column

def weight_function(p, lambda_, T):
    """ Computes w_n(p) for a given m/T = lambda """
    E = np.sqrt(p**2 + (lambda_ * T)**2)  # m = lambda * T
    return np.exp(p / T) / (np.exp(E / T) - 1)

def max_weight_condition(lambda_, p_grid, T):
    """ Finds the maximum value of w_n(p) on the Gauss-Laguerre grid """
    w_values = np.array([weight_function(p, lambda_, T) for p in p_grid])
    return np.max(w_values)

def rational_polynomial_fit(x):
    """ Computes the rational polynomial fit from pion_thermal_weight_max function """
    x2 = x * x
    x3 = x2 * x
    x4 = x3 * x
    numerator = (143206.88623164667 - 95956.76008684626*x - 21341.937407169076*x2 +
                 14388.446116867359*x3 - 6083.775788504437*x4)
    denominator = (-0.3541350577684533 + 143218.69233952634*x - 24516.803600065778*x2 -
                   115811.59391199696*x3 + 35814.36403387459*x4)
    return 1.00001 * (numerator / denominator)

def fit_rational_polynomial(lambda_values, max_w_values):
    """ Finds the coefficients for a rational polynomial fit of max w_n(p) """
    from scipy.optimize import curve_fit
    
    def rational_poly(x, a0, a1, a2, a3, a4, b0, b1, b2, b3, b4):
        numerator = a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4
        denominator = b0 + b1*x + b2*x**2 + b3*x**3 + b4*x**4
        return numerator / denominator
    
    popt, _ = curve_fit(rational_poly, lambda_values, max_w_values, maxfev=10000)
    return popt, rational_poly

def find_max_as_function_of_lambda(filepath, T=1):
    """ Finds max w_n(p) as a function of m/T and fits a polynomial """
    p_grid = read_momentum_grid(filepath)
    lambda_values = np.linspace(0.1, 0.8554, 50)
    max_w_values = [max_weight_condition(lambda_, p_grid, T) for lambda_ in lambda_values]
    
    # Fit rational polynomial to computed max_w_values
    coeffs, rational_poly = fit_rational_polynomial(lambda_values, max_w_values)
    fit_values = rational_poly(lambda_values, *coeffs)
    ref_values = [rational_polynomial_fit(lmb) for lmb in lambda_values]
    
    print("Fitted rational polynomial coefficients:", coeffs)
    
    # Triple plot comparison
    plt.figure(figsize=(8,6))
    plt.scatter(lambda_values, max_w_values, label="Computed max w_n(p)", color='blue')
    plt.plot(lambda_values, ref_values, label="Reference fit", linestyle='--', color='red')
    plt.plot(lambda_values, fit_values, label="Fitted polynomial", linestyle='-.', color='green')
    plt.axhline(y=1, color='black', linestyle='--')
    plt.xlabel("m/T")
    plt.ylabel("Max w_n(p)")
    plt.legend()
    plt.title("Comparison of Computed, Reference, and Fitted max w_n(p)")
    plt.grid()
    plt.savefig("max_w_comparison.png")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_gauss_laguerre_grid>")
        sys.exit(1)
    
    filepath = sys.argv[1]  # Get file path from command line argument
    find_max_as_function_of_lambda(filepath)
