import phoebe
from phoebe import u
import numpy as np

# Main execution block for safe multiprocessing
# The if __name__ == '__main__': guard ensures that the script's main execution is only run when the script is executed directly,
# preventing child processes from re-importing the module during multiprocessing, which is critical for the emcee solver.
# This is particularly important on Windows, where the 'spawn' start method is used for multiprocessing.
# Expected result: Safely executes the enclosed code, allowing emcee to run without errors. No return value.
if __name__ == '__main__':
    # Create bundle for detached binary
    # phoebe.default_binary(): Creates a default binary system with two stars and an orbit. Type: phoebe.frontend.bundle.Component (binary object).
    # The call builds a hierarchical bundle structure for a binary star system with dynamical properties.
    # Expected result: A PHOEBE Bundle object 'b' representing the system, ready for parameter setting and computations. Type of b: phoebe.frontend.bundle.Bundle.
    b = phoebe.default_binary()

    # Set parameters
    # b['period@orbit'] = 1.0 * u.day: Sets the orbital period parameter in the bundle. Parameters passed: key='period@orbit' (str), value=1.0 * u.day (astropy.units.Quantity with unit day).
    # Expected result: Updates the 'period' parameter in the orbit component to 1.0 days. No return value (in-place modification).
    b['period@orbit'] = 1.0 * u.day

    # b['incl@orbit'] = 85.0 * u.deg: Sets the inclination angle. Parameters passed: key='incl@orbit' (str), value=85.0 * u.deg (astropy.units.Quantity with unit degree).
    # Expected result: Updates the 'incl' parameter in the orbit component to 85 degrees. No return value.
    b['incl@orbit'] = 85.0 * u.deg

    # b['teff@primary'] = 6000.0 * u.K: Sets the effective temperature for the primary star. Parameters passed: key='teff@primary' (str), value=6000.0 * u.K (astropy.units.Quantity with unit Kelvin).
    # Expected result: Updates 'teff' in the primary star component. No return value.
    b['teff@primary'] = 6000.0 * u.K

    # b['teff@secondary'] = 5000.0 * u.K: Sets the effective temperature for the secondary star. Parameters passed: key='teff@secondary' (str), value=5000.0 * u.K (astropy.units.Quantity with unit Kelvin).
    # Expected result: Updates 'teff' in the secondary star component. No return value.
    b['teff@secondary'] = 5000.0 * u.K

    # b['requiv@primary'] = 1.5 * u.solRad: Sets the equivalent radius for the primary star. Parameters passed: key='requiv@primary' (str), value=1.5 * u.solRad (astropy.units.Quantity with unit solar radius).
    # Expected result: Updates 'requiv' in the primary star component. No return value.
    b['requiv@primary'] = 1.5 * u.solRad

    # b['requiv@secondary'] = 1.2 * u.solRad: Sets the equivalent radius for the secondary star. Parameters passed: key='requiv@secondary' (str), value=1.2 * u.solRad (astropy.units.Quantity with unit solar radius).
    # Expected result: Updates 'requiv' in the secondary star component. No return value.
    b['requiv@secondary'] = 1.2 * u.solRad

    # b['ecc@orbit'] = 0.0: Sets the eccentricity of the orbit. Parameters passed: key='ecc@orbit' (str), value=0.0 (float).
    # Expected result: Updates 'ecc' in the orbit component to 0 (circular orbit). No return value.
    b['ecc@orbit'] = 0.0

    # b['irrad_frac_refl_bol@primary'] = 0.0: Sets the irradiation reflection fraction for the primary star. Parameters passed: key='irrad_frac_refl_bol@primary' (str), value=0.0 (float).
    # Expected result: Disables irradiation reflection for primary. No return value.
    b['irrad_frac_refl_bol@primary'] = 0.0

    # b['irrad_frac_refl_bol@secondary'] = 0.0: Sets the irradiation reflection fraction for the secondary star. Parameters passed: key='irrad_frac_refl_bol@secondary' (str), value=0.0 (float).
    # Expected result: Disables irradiation reflection for secondary. No return value.
    b['irrad_frac_refl_bol@secondary'] = 0.0

    # b['atm@primary'] = 'blackbody': Sets the atmosphere model for the primary star. Parameters passed: key='atm@primary' (str), value='blackbody' (str).
    # Expected result: Applies blackbody atmosphere to the primary star. No return value.
    b['atm@primary'] = 'blackbody'

    # b['atm@secondary'] = 'blackbody': Sets the atmosphere model for the secondary star. Parameters passed: key='atm@secondary' (str), value='blackbody' (str).
    # Expected result: Applies blackbody atmosphere to the secondary star. No return value.
    b['atm@secondary'] = 'blackbody'

    # Add synthetic LC dataset (simulate times and fluxes with noise)
    times = np.linspace(0, 1, 101)

    # b.add_dataset('lc', times=times, passband='Johnson:V', dataset='lc01'): Adds a light curve dataset to the bundle.
    # Parameters passed: kind='lc' (str), times=times (numpy.ndarray of float, length 101), passband='Johnson:V' (str), dataset='lc01' (str).
    # Expected result: Creates a new dataset component 'lc01' in the bundle for light curve computations. Returns the dataset twig (str).
    b.add_dataset('lc', times=times, passband='Johnson:V', dataset='lc01')

    # b['ld_mode@lc01@primary'] = 'manual': Sets the limb darkening mode for the primary star in the light curve dataset. Parameters passed: key='ld_mode@lc01@primary' (str), value='manual' (str).
    # Expected result: Enables manual specification of limb darkening coefficients for the primary star in the light curve. No return value.
    b['ld_mode@lc01@primary'] = 'manual'

    # b['ld_mode@lc01@secondary'] = 'manual': Sets the limb darkening mode for the secondary star in the light curve dataset. Parameters passed: key='ld_mode@lc01@secondary' (str), value='manual' (str).
    # Expected result: Enables manual specification of limb darkening coefficients for the secondary star in the light curve. No return value.
    b['ld_mode@lc01@secondary'] = 'manual'

    # b['ld_func@lc01@primary'] = 'linear': Sets the limb darkening functional form for the primary star in the light curve. Parameters passed: key='ld_func@lc01@primary' (str), value='linear' (str).
    # Expected result: Applies linear limb darkening law to the primary star in the light curve dataset. No return value.
    b['ld_func@lc01@primary'] = 'linear'

    # b['ld_func@lc01@secondary'] = 'linear': Sets the limb darkening functional form for the secondary star in the light curve. Parameters passed: key='ld_func@lc01@secondary' (str), value='linear' (str).
    # Expected result: Applies linear limb darkening law to the secondary star in the light curve dataset. No return value.
    b['ld_func@lc01@secondary'] = 'linear'

    # b['ld_coeffs@lc01@primary'] = [0.5]: Sets limb darkening coefficients for the primary star in the light curve dataset. Parameters passed: key='ld_coeffs@lc01@primary' (str), value=[0.5] (list of float).
    # Expected result: Updates limb darkening coefficients for the primary star in 'lc01'. No return value.
    b['ld_coeffs@lc01@primary'] = [0.5]

    # b['ld_coeffs@lc01@secondary'] = [0.5]: Sets limb darkening coefficients for the secondary star in the light curve dataset. Parameters passed: key='ld_coeffs@lc01@secondary' (str), value=[0.5] (list of float).
    # Expected result: Updates limb darkening coefficients for the secondary star in 'lc01'. No return value.
    b['ld_coeffs@lc01@secondary'] = [0.5]

    # b['gravb_bol@primary'] = 0.32: Sets the bolometric gravity brightening coefficient for the primary star. Parameters passed: key='gravb_bol@primary' (str), value=0.32 (float).
    # Expected result: Updates gravity brightening for the primary star to value appropriate for convective atmosphere. No return value.
    b['gravb_bol@primary'] = 0.32

    # b['gravb_bol@secondary'] = 0.32: Sets the bolometric gravity brightening coefficient for the secondary star. Parameters passed: key='gravb_bol@secondary' (str), value=0.32 (float).
    # Expected result: Updates gravity brightening for the secondary star to value appropriate for convective atmosphere. No return value.
    b['gravb_bol@secondary'] = 0.32

    # Compute model fluxes (true model)
    # b.run_compute(model='reference_model'): Runs the forward model computation.
    # Parameters passed: model='reference_model' (str, labels the model).
    # Expected result: Computes synthetic fluxes and stores them in the bundle under 'reference_model'. Returns a ParameterSet of computed values (phoebe.parameters.ParameterSet).
    b.run_compute(model='reference_model')

    # Add Gaussian noise to simulate data
    # b['fluxes@reference_model@lc01'].value: Accesses the computed fluxes from the model.
    # Parameters passed: key='fluxes@reference_model@lc01' (str, twig specifying component, model, and dataset).
    # Expected result: Returns the flux array as a numpy.ndarray of float (length 101). The .value attribute extracts the numerical values.
    fluxes = b['fluxes@reference_model@lc01'].value + 0.01 * np.random.randn(len(times))

    # b['fluxes@lc01@dataset'] = fluxes: Sets the observed fluxes in the dataset.
    # Parameters passed: key='fluxes@lc01@dataset' (str), value=fluxes (numpy.ndarray of float).
    # Expected result: Updates the 'fluxes' parameter in the 'lc01' dataset component. No return value.
    b['fluxes@lc01@dataset'] = fluxes

    # b['sigmas@lc01@dataset'] = 0.01 * np.ones_like(fluxes): Sets the flux uncertainties.
    # Parameters passed: key='sigmas@lc01@dataset' (str), value=0.01 * np.ones_like(fluxes) (numpy.ndarray of float, all 0.01).
    # Expected result: Updates the 'sigmas' parameter in the 'lc01' dataset component. No return value.
    b['sigmas@lc01@dataset'] = 0.01 * np.ones_like(fluxes)

    # Set priors on the parameters we want to fit.
    # To show covariance, we will fit for both teffs.
    # We will create one distribution collection to hold all our priors.
    b.add_distribution('teff@primary', phoebe.gaussian(6000, 500), distribution='mcmc_priors')
    b.add_distribution('teff@secondary', phoebe.gaussian(5000, 500), distribution='mcmc_priors')

    # Add emcee solver
    # b.add_solver('emcee', solver='emcee_solver'): Adds an emcee MCMC solver to the bundle.
    # Parameters passed: kind='emcee' (str), solver='emcee_solver' (str, labels the solver).
    # Expected result: Creates a solver component 'emcee_solver' in the bundle. Returns the solver twig (str).
    b.add_solver('emcee', solver='emcee_solver')

    # Set initialization for emcee to use our prior distribution collection.
    b['init_from@emcee_solver'] = 'mcmc_priors'

    # Run emcee sampler (increased iterations for convergence)
    # The 'priors' argument tells the solver which parameters to fit.
    emcee_sol = b.run_solver('emcee_solver',
                             solution='emcee_sol',
                             nwalkers=32,
                             niters=500,
                             priors='mcmc_priors',
                             compute='phoebe01',
                             dataset='lc01')

    # Check convergence (trace plot; in full run, inspect for stationarity)
    # This plot is essential for diagnosing whether the MCMC has converged.
    b.plot(solution='emcee_sol', style='trace', show=True)

    # Adopt posterior distributions (without adopting values)
    # b.adopt_solution('emcee_sol', adopt_values=False, adopt_distributions=True, distribution='myposterior', burnin=100, thin=5, lnprob_cutoff=None): Adopts posteriors from sampler.
    # Parameters passed: solution='emcee_sol' (str), adopt_values=False (bool), adopt_distributions=True (bool), distribution='myposterior' (str), burnin=100 (int, discards initial iterations), thin=5 (int, samples every 5th iteration), lnprob_cutoff=None (NoneType).
    # Expected result: Attaches posterior distributions to parameters under collection 'myposterior'. Returns a list of adopted parameters (list of str).
    b.adopt_solution('emcee_sol',
                     adopt_values=False,
                     adopt_distributions=True,
                     distribution='myposterior',
                     burnin=100,
                     thin=5,
                     lnprob_cutoff=None)

    # Access and print posterior for teff@primary
    param_pri = b.get_parameter(qualifier='teff', component='primary', context='component')
    posterior_dist_pri = param_pri.get_distribution(distribution='myposterior')
    print(f"Posterior median for teff@primary: {posterior_dist_pri.median():.2f} K")
    sigma_pri = (posterior_dist_pri.uncertainties()[2] - posterior_dist_pri.uncertainties()[0]) / 2
    print(f"Posterior 1-sigma uncertainty for teff@primary: {sigma_pri:.2f} K")

    # Access and print posterior for teff@secondary
    param_sec = b.get_parameter(qualifier='teff', component='secondary', context='component')
    posterior_dist_sec = param_sec.get_distribution(distribution='myposterior')
    print(f"Posterior median for teff@secondary: {posterior_dist_sec.median():.2f} K")
    sigma_sec = (posterior_dist_sec.uncertainties()[2] - posterior_dist_sec.uncertainties()[0]) / 2
    print(f"Posterior 1-sigma uncertainty for teff@secondary: {sigma_sec:.2f} K")


    # For joint posteriors, access collection
    # b.get_distribution_collection(distribution='myposterior'): Retrieves a collection of all posteriors.
    # Parameters passed: distribution='myposterior' (str).
    # Expected result: Returns a DistributionCollection object. Type of dist_coll: phoebe.dependencies.distl.DistributionCollection.
    dist_coll = b.get_distribution_collection(distribution='myposterior')

    # The corner plot is the best way to visualize the joint posterior distributions
    # and see the covariance between the fitted parameters.
    b.plot(distribution='myposterior', style='corner', show=True)

