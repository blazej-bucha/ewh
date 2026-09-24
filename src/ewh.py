# Module to compute surface density (routine "surface_density") and equivalent
# water heights (routine "equivalent_water_height") at the Earth's surface.


import sys
import numpy as np
import scipy as sp
import pyharm as ph


# Enable polar optimization in PyHarm to accelerate spherical-harmonic
# transforms.  The performance gain generally improves with truncation degree
# of spherical-harmonic transforms.
ph.glob.Constants().polar_optimization_a2 = 0.01


def surface_density(shcs_grav,
                    nmax_grav,
                    shcs_shape,
                    nmax_shape,
                    pmax,
                    nmax_density,
                    iterations,
                    M):
    """
    Returns spherical-harmonic coefficients of a surface density function
    residing on the topography from coefficients of the gravitational potential
    and of the Earth's shape.

    Parameters
    ----------
    shcs_grav : A `pyharm.shc.Shc` class instance
        Spherical-harmonic coefficients of the gravitational potential
    nmax_grav : integer
        Degree to truncate `shcs_grav` in the gravity inversion
    shcs_shape : A `pyharm.shc.Shc` class instance
        Spherical-harmonic coefficients of the Earth's shape
    nmax_shape : integer
        Degree to truncate `shcs_shape` in the gravity inversion
    pmax : integer
        Truncation order of a binomial expansion.  The reference radius of the
        binomial expansion is chosen to be the mean sphere, that is,
        `shcs_shape.get_coeffs(n=0, m=0)[0]`
    nmax_density : integer
        Truncation degree of the output spherical-harmonic coefficients of the
        surface density
    iterations : integer
        Number of iterations
    M : floating point number
        Mass of the Earth to scale the spherical-harmonic coefficients of the
        surface density

    Returns
    -------
    out : A `pyharm.shc.Shc` class instance
        Spherical-harmonic coefficients of the surface density
    """


    # Sanity checks
    # -------------------------------------------------------------------------
    if not isinstance(shcs_grav, ph.shc.Shc):
        sys.exit('\'shcs_grav\' must be an instance of \'pyharm.shc.Shc\'.')

    if not isinstance(nmax_grav, int):
        sys.exit('\'nmax_grav\' must be an integer.')

    if nmax_grav < 0:
        sys.exit('\'nmax_grav\' cannot be negative.')

    if not isinstance(shcs_shape, ph.shc.Shc):
        sys.exit('\'shcs_shape\' must be an instance of \'pyharm.shc.Shc\'.')

    if not isinstance(nmax_shape, int):
        sys.exit('\'nmax_shape\' must be an integer.')

    if nmax_shape < 0:
        sys.exit('\'nmax_shape\' cannot be negative.')

    if not isinstance(pmax, int):
        sys.exit('\'pmax\' must be an integer.')

    if pmax < 0:
        sys.exit('\'pmax\' cannot be negative.')

    if not isinstance(iterations, int):
        sys.exit('\'iterations\' must be an integer.')

    if iterations < 0:
        sys.exit('\'iterations\' cannot be negative.')

    if not np.isreal(M):
        sys.exit('\'M\' must be a real number.')
    # -------------------------------------------------------------------------


    # Topographic height function
    # -------------------------------------------------------------------------
    print('Computing spherical-harmonic coefficients of the topographic '
          'height function')

    # Create a copy of shape coefficients that will hold the coefficients of
    # the topographic height function at the end of this code block
    shcs_topo_height = ph.shc.Shc.from_copy(shcs_shape)

    # Reference the topography to the mean sphere
    R = shcs_topo_height.get_coeffs(n=0)[0]
    shcs_topo_height.set_coeffs(n=0, c=np.zeros((1,)))

    # Compute the coefficients of the topographic height function
    shcs_topo_height.c[:] /= R
    shcs_topo_height.s[:] /= R
    # -------------------------------------------------------------------------


    # Initialize spherical-harmonic coefficients of "mu0 = sigma * J" using
    # spherical approximation
    # -------------------------------------------------------------------------
    print('Initializing spherical-harmonic coefficients of the surface '
          'density')

    shcs_mu0_init = ph.shc.Shc.from_copy(shcs_grav,
                                         nmax=min(shcs_grav.nmax,
                                                  nmax_density),
                                         nmax_shcs_out=nmax_density)
    shcs_mu0_init.mu = shcs_mu0_init.r = 1.0

    c = M / (4.0 * np.pi * R**2)
    t = shcs_grav.r / R
    n = np.arange(nmax_density + 1, dtype=np.float64)
    shcs_mu0_init.mul_degree_wise(c * ((2.0 * n + 1.0)) * t**n)
    # -------------------------------------------------------------------------


    # Iterative scheme for "mu0"
    # -------------------------------------------------------------------------
    print('Starting iterations')

    # Initialize "shcs_mu0_prev" to "shcs_mu0_init" from the spherical
    # model
    shcs_mu0_prev = ph.shc.Shc.from_copy(shcs_mu0_init)

    # Create an array of spherical harmonic degrees up to the maximum degree
    # that will be used for grid sampling in spherical-harmonic transforms
    n = np.arange(pmax * nmax_shape + nmax_density + 1, dtype=np.float64)

    print('\tRMS differences between two consecutive iterations:')

    # Iterations
    for i in range(iterations + 1):

        # Compute the sum over "p = 1, .., pmax"
        shcs_mu_sum = ph.shc.Shc.from_zeros(nmax_density)
        for p in range(1, pmax + 1):
            # Maximum degree for the Gauss--Legendre grid for a given "p"
            N_tmp = p * nmax_shape + nmax_density

            # Create Gauss--Legendre grid
            grd = ph.crd.PointGridGL(N_tmp)

            # Synthesize "mu0" using the coefficients from the previous
            # iteration
            mu0 = ph.shs.point(grd, shcs_mu0_prev, nmax_density)

            # Compute the topographic height function "h = (rs - R) / R", where
            # "rs" is the full spherical radius to the Earth's shape
            h = ph.shs.point(grd, shcs_topo_height, nmax_shape)

            # Compute spherical-harmonic coefficients of "h**p * mu0"
            shcs_mu_p = ph.sha.point(grd, h**p * mu0, N_tmp)

            # Apply the "n + 2 choose p" factor
            cnp = sp.special.binom(n[:(N_tmp + 1)] + 2.0, p)
            shcs_mu_p.mul_degree_wise(cnp)

            # Add the contribution of the "p"th power in "shcs_mu_p" to the sum
            # over "p = 1, ..., pmax" in "shcs_mu_sum"
            shcs_mu_sum += shcs_mu_p

        # Subtract "shcs_mu_sum" from the coefficients of the initial spherical
        # model
        shcs_mu0  = ph.shc.Shc.from_copy(shcs_mu0_init)
        shcs_mu0 -= shcs_mu_sum

        # Compute RMS between the current and previous approximation of the
        # surface density to reveal convergence rate
        rms = np.sqrt(shcs_mu0.ddv(shcs_mu0_prev).sum())

        print(f'\titeration: {i}, rms: {rms} kg * m**-2')

        # Use the new coefficients of "mu0" as old ones in the next iteration
        # cycle
        shcs_mu0_prev = ph.shc.Shc.from_copy(shcs_mu0)
    # -------------------------------------------------------------------------


    # Get coefficients of "sigma" from "mu0"
    # -------------------------------------------------------------------------
    # Create Gauss--Legendre grid associated with the degree "nmax_density".
    #
    # Degree higher than "nmax_density" could be used here to suppress
    # numerical errors if "nmax_density" is too low
    grd = ph.crd.PointGridGL(nmax_density)

    # Synthesize "mu0 = sigma * J"
    mu0 = ph.shs.point(grd, shcs_mu0, shcs_mu0.nmax)

    # Compute the spherical radius at "grd"
    rs  = ph.shs.point_guru(grd, shcs_shape, nmax_shape, 0, 0, 0)

    # Compute "dr / dlat".  Remember that "grd.r" is "1.0" for all latitudes,
    # so "1 / r * dr / dlat" indeed reduces to "dr / dlat".
    rsp = ph.shs.point_guru(grd, shcs_shape, nmax_shape, 0, 1, 0)

    # Compute "1 / (cos(lat)) * dr / dlon".  The same comment in relation to
    # "grd.r" applies here.
    rsl = ph.shs.point_guru(grd, shcs_shape, nmax_shape, 0, 0, 1)

    # Compute the topography related factor that ensures we get the surface
    # density at the topography
    J   = np.sqrt(1.0 + (rsp**2 + rsl**2) / rs**2)

    # Compute and return spherical-harmonic coefficients of "sigma = mu0 / J"
    return ph.sha.point(grd, mu0 / J, grd.nmax)
    # -------------------------------------------------------------------------


def love_numbers(nmax):
    """
    Returns a numpy array of load Love numbers up to degree `nmax`.

    Parameters
    ----------
    nmax : integer
        Maximum harmonic degree

    Returns
    -------
    out : numpy array
        Load Love numbers up to degree `nmax`
    """

    # Love coefficients from Wahr and Molenaar (1998)
    n = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 40, 50, \
                  70, 100, 150, 200])
    kn = np.array([0.0, 0.027, -0.303, -0.194, -0.132, -0.104, -0.089, \
                   -0.081, -0.076, -0.072, -0.069, -0.064, -0.058, -0.051, \
                   -0.040, -0.033, -0.027, -0.020, -0.014, -0.010, -0.007])

    # Piecewise cubic interpolator
    interp = sp.interpolate.CubicSpline(n, kn)

    return interp(np.arange(nmax + 1))


def equivalent_water_height(sigma):
    """
    Transforms spherical-harmonic coefficients of a surface density returned by
    the `surface_density` routine to coefficients of equivalent water heights.

    Parameters
    ----------
    sigma : A `pyharm.shc.Shc` class instance
        Spherical-harmonic coefficients of the surface density.  These
        coefficients will be overwritten by the coefficients of equivalent
        water heights.

    Returns
    -------
    None
    """

    sigma.mul_degree_wise(1.0 / ((1 + love_numbers(sigma.nmax)) * 1000.0))

    return

