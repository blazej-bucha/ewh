# This script shows how to determine equivalent water heights at the Earth's
# surface from a single GRACE gravitational model.


# IMPORTS
# =============================================================================
# Import PyHarm, a C/Python library for spherical-harmonic transforms (see
# "https://charmlib.org")
import pyharm as ph

# Import "ewh.py", which contains routines to compute the equivalent water
# heights
import ewh

# Import numpy for array calculations
import numpy as np

# Import matplotlib for plotting
import matplotlib.pyplot as plt
# =============================================================================


# INPUTS
# =============================================================================
# Path to spherical-harmonic coefficients of the gravitational potential from
# GRACE
SHCS_GRAV = '../data/input/GSM-2_2002095-2002120_GRAC_UTCSR_BA01_0600.gfc'

# Truncation degree of the GRACE model
NMAX_GRAV = 60

# Path to spherical-harmonic coefficients of the Earth's shape
SHCS_SHAPE = '../data/input/earth-shape.tbl'

# Truncation degree of the Earth shape model
NMAX_SHAPE = 60

# Maximum topography power for a binomial expansion
MAX_SHAPE_POWER = 10

# Equivalent water heights will be recovered up to degree "NMAX_GRAV
# + EXTRA_SHCS".
#
# "EXTRA_SHCS" should be large enough to capture signals in equivalent water
# heights that are found beyond the truncation degree of the gravitational
# potential for non-spherical Earth shape models.
EXTRA_SHCS = 264

# Number of iterations
ITERATIONS = 10

# Newton's gravitational constant ("m**3 * kg**-1 * s**-2")
G = 6.67430 * 10**-11

# Name of the output file with the spherical-harmonic coefficients of
# equivalent water heights
OUTPUT_FILE_SHCS = '../data/output/ewh.tbl'

# Name of the output file showing the equivalent water heights in the spatial
# domain
OUTPUT_FILE_FIGURE = '../data/output/ewh.png'
# =============================================================================


# SPHERICAL-HARMONIC COEFFICIENTS OF EQUIVALENT WATER HEIGHTS
# =============================================================================
# Read GRACE spherical-harmonic coefficients
print(f'Reading potential coefficients from \'{SHCS_GRAV}\'')
shcs_grav = ph.shc.Shc.from_file('gfc', SHCS_GRAV, NMAX_GRAV)

# Read Earth shape spherical-harmonic coefficients
print(f'Reading Earth shape coefficients from \'{SHCS_SHAPE}\'')
shcs_shape = ph.shc.Shc.from_file('tbl', SHCS_SHAPE, NMAX_SHAPE)

# Mass of the Earth
M = shcs_grav.mu / G

# Truncation degree of the surface density and of equivalent water heights
nmax_density = shcs_grav.nmax + EXTRA_SHCS

# Compute spherical-harmonic coefficients of the surface density.  Later, these
# coefficients are transformed into the coefficients of equivalent water
# heights, so named them "shcs_ewh".
print('Computing surface density')
shcs_ewh = ewh.surface_density(shcs_grav,
                               NMAX_GRAV,
                               shcs_shape,
                               NMAX_SHAPE,
                               MAX_SHAPE_POWER,
                               nmax_density,
                               ITERATIONS,
                               M)

# Transform the coefficients of the surface density to coefficients of
# equivalent water heights
print('Computing equivalent water heights')
ewh.equivalent_water_height(shcs_ewh)

# Write spherical-harmonic coefficients of the equivalent water heights to
# a text file
print(f'Writing coefficients of equivalent water heights to '
      f'\'{OUTPUT_FILE_SHCS}\'')
shcs_ewh.to_file('tbl', OUTPUT_FILE_SHCS)
# =============================================================================


# EQUIVALENT WATER HEIGHTS IN THE SPATIAL DOMAIN
# =============================================================================
print('Computing equivalent water heights at global grid')

# For simplicity, the equivalent water heights will be computed at
# a Driscoll--Healy point grid
grd = ph.crd.PointGridDH2(nmax_density)

# Synthesize of the equivalent water heights
ewh_grd = ph.shs.point(grd, shcs_ewh, shcs_ewh.nmax)

# Simple plotting
print('Plotting equivalent water heights')
latmin = np.rad2deg(grd.lat.min())
latmax = np.rad2deg(grd.lat.max())
lonmin = np.rad2deg(grd.lon.min())
lonmax = np.rad2deg(grd.lon.max())
plt.imshow(ewh_grd,
           extent=(lonmin, lonmax, latmin, latmax))
plt.xlabel('Longitude (deg)')
plt.ylabel('Latitude (deg)')
plt.colorbar()
plt.title('Equivalent water heights (m)')
plt.savefig(OUTPUT_FILE_FIGURE, dpi=600)
plt.close('all')
# =============================================================================


# =============================================================================
print('\nDone!')
print(f'Spherical-harmonic coefficients of equivalent water heights are '
      f'stored in \'{OUTPUT_FILE_SHCS}\'.')
print(f'The equivalent water heights are plotted in '
      f'\'{OUTPUT_FILE_FIGURE}\'.  Keep in mind that the figure shows a '
      f'snapshot for a single month.  To get, for instance, linear trends, '
      f'you\'d have to analyze snapshots like this one from every monthly '
      f'solution as a time-series.')
# =============================================================================
