"""
define some constants and conversions

"""

from astropy.cosmology import Planck15 as plk
from astropy.constants import c, k_B, m_p, m_e, h  # pyright: ignore[reportAttributeAccessIssue]
import numpy as np
import astropy.units as u

# physical constants
PI = np.pi
C_PLK = h  # Planck constant
C_BOL = k_B  # Boltzmann constant

# cosmology constants
LIGHT_SPEED = c.to(u.Unit("km/s"))  # speed of light in the vaccumn
C_HUB = (
    plk.h  # pyright: ignore[reportAttributeAccessIssue]
)  # dimensionless hubble constant

# HI constants
A_HI = u.Quantity(2.86888e-15, "1/s")  # HI spontaneous emission rate
HI_MASS = m_p + m_e  # Atomatic HI mass
HI_REST_FREQ = u.Quantity(
    1420405751.78, "Hz"
)  # HI 21cm emission frequency in the rest frame
HI_REST_WAVE = HI_REST_FREQ.to(
    u.Unit("cm"), equivalencies=u.spectral()
)  # HI 21cm emission wavelength in the rest frame

# conversion constants
SIGMA2FWHM = (
    8 * np.log(2)
) ** 0.5  # turn guass sigma into full width at half maximum power point
FWHM2SIGMA = 1 / SIGMA2FWHM
