import numpy as np
from D2Q9 import D2Q9  # Assuming the D2Q9 class is defined in a file named D2Q9.py


def test_main():
    """Test the main function of D2Q9.py using assert statements."""

    # Expected input data from the D2Q9 class
    nodes = {
        "A": np.array([1.63, 0.61, 0.41, 0.27, 0.41, 0.15, 0.07, 0.07, 0.16]),
        "B": np.array([1.67, 0.42, 0.42, 0.42, 0.42, 0.1, 0.11, 0.1, 0.11]),
        "C": np.array([1.66, 0.5, 0.42, 0.35, 0.42, 0.12, 0.09, 0.08, 0.13]),
    }

    # Initialize the D2Q9 lattice
    d2q9 = D2Q9()

    # Expected results (replace with actual expected values based on calculations or specifications)
    expected_results = {
        "A": {
            "density": 3.78,
            "velocity": np.array([0.1349, -0.0026]),
            "equilibrium": np.array([1.634109, 0.61293, 0.405207, 0.272932, 
                                    0.411874, 0.152066, 0.067740, 0.068732, 0.154407]),
        },

        "B": {
            "density": 3.77,
            "velocity": np.array([0.0, 0.3681e-17]),
            "equilibrium": np.array([1.675556, 0.418889, 0.418889, 0.418889,
                                    0.418889, 0.104722, 0.104722, 0.104722, 0.104722]),
        },

        "C": {
            "density": 3.77,
            "velocity": np.array([0.0610, 0.0]),
            "equilibrium": np.array([1.666201, 0.500233, 0.416550, 0.346899,
                                    0.416550, 0.125058, 0.086725, 0.086725, 0.125058]),
        },
    }

    # Test each node
    for node, f in nodes.items():
        rho, u = d2q9.moment_rho_u(f)
        feq = d2q9.compute_equilibrium(rho, u)

        # Validate density
        assert np.isclose(rho, expected_results[node]["density"], atol=1e-2), f"Density mismatch for node {node}"  # Tolerance level

        # Validate velocity
        assert np.allclose(u, expected_results[node]["velocity"], atol=1e-2), f"Velocity mismatch for node {node}"

        # Validate equilibrium distribution
        assert np.allclose(feq, expected_results[node]["equilibrium"], atol=1e-3), f"Equilibrium distribution mismatch for node {node}"

    print("All tests passed successfully!")


if __name__ == "__main__":
    test_main()
