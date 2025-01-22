import numpy as np
from D2Q9 import D2Q9

def test_equilibrium_distribution():
    """
    Test the equilibrium distribution function with varying inputs for accuracy.
    """
    # Instantiate the D2Q9Lattice class
    lattice = D2Q9()

    # Define test cases with inputs and expected outputs
    test_cases = [
        {
            "rho": 1.0,
            "ux": 0.0,
            "uy": 0.0,
            "expected": np.round(np.array([
                4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36
            ]), 3)
        },
        {
            "rho": 1.2,
            "ux": 0.1,
            "uy": 0.0,
            "expected": np.round(np.array([
                0.525, 0.177, 0.131, 0.097, 0.131, 0.044, 0.024, 0.024, 0.044
            ]), 3)
        },
        {
            "rho": 0.8,
            "ux": 0.0,
            "uy": -0.1,
            "expected": np.round(np.array([
                0.35, 0.088, 0.065, 0.088, 0.118, 0.016, 0.016, 0.03, 0.03
            ]), 3)
        }
    ]

    # Run each test case
    for i, case in enumerate(test_cases):
        rho = case["rho"]
        ux = case["ux"]
        uy = case["uy"]
        expected = case["expected"]

        # Compute equilibrium distribution using the D2Q9Lattice class method
        result = lattice.equilibrium_distribution_single_unit(rho, ux, uy)

        # Compare result to expected value
        if np.allclose(result, expected):
            print(f"Test case {i+1} passed.")
        else:
            print(f"Test case {i+1} failed.\nExpected: {expected}\nGot: {result}")

def test_initialize_single_lattice_unit():
    """
    Test the initialization of a single lattice unit.
    """
    # Instantiate the D2Q9Lattice class
    lattice = D2Q9()

    # Expected output for uniform density (rho = 1.0)
    expected = np.round(np.array([
        4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36
    ]), 3)

    # Initialize lattice unit using the D2Q9Lattice class method
    result = lattice.initialize_single_lattice_unit()[0]  # Extract the single lattice unit values

    # Compare result to expected value
    if np.allclose(result, expected):
        print("Initialization test passed.")
    else:
        print(f"Initialization test failed.\nExpected: {expected}\nGot: {result}")

# Run the tests
if __name__ == "__main__":
    print("Running tests for the D2Q9 lattice...")
    test_initialize_single_lattice_unit()
    test_equilibrium_distribution()
