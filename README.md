# PHOEBE 2 Example: MCMC Fitting with emcee

This script provides a comprehensive, educational example of how to use PHOEBE 2. It demonstrates how to:
1.  Model a detached binary star system.
2.  Generate a synthetic light curve with simulated noise.
3.  Perform parameter estimation using the `emcee` Markov Chain Monte Carlo (MCMC) solver.
4.  Analyze the results to determine the posterior probability of a model parameter.

## Prerequisites

To run this script, you will need Python with the following libraries installed:
*   `phoebe`
*   `numpy`

You can typically install these using pip:
```bash
pip install phoebe numpy
```

## How to Run

To execute the script, simply run the following command from your terminal in the `example` directory:

```bash
python example.py
```

---

## A Note on Multiprocessing (`if __name__ == '__main__':`)

You will notice the entire script is wrapped in an `if __name__ == '__main__':` block. This is a critical feature for ensuring the code runs correctly, especially on Windows. The `emcee` solver uses multiprocessing to run its calculations in parallel. This guard prevents child processes from re-importing and re-executing the script, which would lead to errors and infinite loops. It is a standard best practice for Python scripts that use multiprocessing.

---

## Workflow Overview

The script follows these key steps:

1.  **System Setup**: A default binary star system is created using `phoebe.default_binary()`. The script then configures the system by setting various physical and computational parameters for the orbit and the two stars.

    The following parameters are set:
    *   **Orbital Parameters**:
        *   `period@orbit`: 1.0 day (Orbital period)
        *   `incl@orbit`: 85.0 degrees (Inclination of the orbital plane)
        *   `ecc@orbit`: 0.0 (Eccentricity, for a circular orbit)
    *   **Stellar Parameters**:
        *   `teff@primary`: 6000 K (Effective temperature of the primary star)
        *   `teff@secondary`: 5000 K (Effective temperature of the secondary star)
        *   `requiv@primary`: 1.5 Solar radii (Equivalent radius of the primary star)
        *   `requiv@secondary`: 1.2 Solar radii (Equivalent radius of the secondary star)
    *   **Atmosphere & Irradiation**:
        *   `atm@primary` / `atm@secondary`: 'blackbody' (Sets a simple blackbody atmospheric model)
        *   `irrad_frac_refl_bol@primary` / `irrad_frac_refl_bol@secondary`: 0.0 (Disables reflection due to irradiation)
        *   `gravb_bol@primary` / `gravb_bol@secondary`: 0.32 (Bolometric gravity brightening coefficient)
    *   **Limb Darkening (for the 'lc01' dataset)**:
        *   `ld_mode`: 'manual' (Enables manual setting of coefficients)
        *   `ld_func`: 'linear' (Uses a linear limb-darkening law)
        *   `ld_coeffs`: [0.5] (The linear limb-darkening coefficient)

2.  **Synthetic Data Generation**:
    *   A synthetic light curve (`lc`) dataset is added to the model with a specified set of observation times.
    *   The forward model is computed using `b.run_compute()` to generate the "true" fluxes for the defined system.
    *   Gaussian noise is added to these fluxes to simulate a realistic observational dataset. The uncertainties (`sigmas`) are also defined.

3.  **Prior Definition**:
    *   A **prior** distribution is defined for the parameter we want to fit. A prior represents our initial belief about a parameter's value before we consider the data. In this example, a Gaussian prior is placed on `teff@primary`.

4.  **MCMC Sampling with `emcee`**:
    *   The `emcee` solver is added to the PHOEBE bundle.
    *   The solver is configured to initialize its "walkers" from the defined prior distribution.
    *   The MCMC simulation is executed with `b.run_solver()`. This explores the parameter space to find the values that best fit the synthetic data.

5.  **Posterior Analysis**:
    *   The results of the MCMC simulation are the **posterior** distributions. A posterior represents our updated belief about a parameter after fitting to the data.
    *   These posteriors are adopted into a new collection using `b.adopt_solution()`. The `burnin` and `thin` parameters are used to discard initial unstable samples and reduce autocorrelation, ensuring the final samples are from the true posterior distribution.
    *   The script then accesses the posterior for `teff@primary` and calculates its median and 1-sigma uncertainty.

## Expected Output

When the script finishes, it will print the statistical properties of the posterior distribution for the primary star's effective temperature. This tells us the most probable value for the parameter and our confidence in that value, given the data. The output will look similar to this:

```
Posterior median for teff@primary: 6001.23 K
Posterior 1-sigma uncertainty: 45.67 K
```

*(Note: The exact numerical values will vary slightly on each run due to the random noise generation and the stochastic nature of the MCMC sampler.)*

## Visualizing the Results

The script includes commented-out plotting commands. It is highly recommended to uncomment these lines to visualize the MCMC results, which is a crucial step in any MCMC analysis.

*   `b.plot(solution='emcee_sol', style='trace', show=True)`: This creates a **trace plot**. It shows the value of the parameter at each step in the chain for each walker. You should look for the traces to be stable and well-mixed, indicating the sampler has converged.
*   `b.plot(distribution='myposterior', style='corner', show=True)`: This creates a **corner plot**. It shows the 1D and 2D posterior distributions for all fitted parameters. This is useful for visualizing the shape of the probability distribution and seeing any correlations between parameters.
