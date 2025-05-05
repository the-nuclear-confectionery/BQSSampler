import numpy as np
import scipy.special as sp
import pandas as pd

def gauss_legendre(n):
    """
    Compute the Gauss-Legendre roots and weights for the given number of points n.
    
    Parameters:
    n (int): Number of quadrature points.
    
    Returns:
    roots (ndarray): The Gauss-Legendre roots.
    weights (ndarray): The Gauss-Legendre weights.
    """
    roots, weights = np.polynomial.legendre.leggauss(n)
    return roots, weights

def gauss_laguerre(a, n):
    """
    Compute the Generalized Gauss-Laguerre roots and weights for the given number of points n and order a.
    
    Parameters:
    a (float): The order parameter alpha.
    n (int): Number of quadrature points.
    
    Returns:
    roots (ndarray): The Gauss-Laguerre roots.
    weights (ndarray): The Gauss-Laguerre weights.
    """
    roots, weights = sp.roots_genlaguerre(n, a)
    return roots, weights

def save_to_dat(filename, roots, weights):
    """
    Save the roots and weights to a .dat file using pandas.
    
    Parameters:
    filename (str): The output file name.
    roots (ndarray): The roots.
    weights (ndarray): The weights.
    """
    # Sort roots and weights based on roots from lowest to highest
    sorted_indices = np.argsort(roots)
    roots_sorted = roots[sorted_indices]
    weights_sorted = weights[sorted_indices]
    
    # Create a DataFrame and save it to .dat
    df = pd.DataFrame({
        "Root": roots_sorted,
        "Weight": weights_sorted
    })
    df.to_csv(filename, sep=' ', index=False, header=False)  # Save with space separator and no header/index
    print(f"Table saved as {filename}")

# Number of points for the quadrature
n = 32  # Example: 5 points for each quadrature

# Gauss-Legendre
roots_legendre, weights_legendre = gauss_legendre(n)
save_to_dat("gauss_legendre.dat", roots_legendre, weights_legendre)

# Gauss-Laguerre with a=0
a_0 = 0
roots_laguerre_0, weights_laguerre_0 = gauss_laguerre(a_0, n)
save_to_dat("gauss_laguerre_a0.dat", roots_laguerre_0, weights_laguerre_0)

# Gauss-Laguerre with a=1
a_1 = 1
roots_laguerre_1, weights_laguerre_1 = gauss_laguerre(a_1, n)
save_to_dat("gauss_laguerre_a1.dat", roots_laguerre_1, weights_laguerre_1)

print("All tables saved as .dat files.")
