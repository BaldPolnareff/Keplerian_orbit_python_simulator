import numpy as np
import scipy as sp

## All the constants are in SI units

# Speed of light
c = 299792458 # m/s

# Earth's gravitational constant 
mu_earth = 3.986004418e14 # m^3/s^2

# Earth's semi-major axis (WGS84)
a_earth = 6378137 # m

# Earth's semi-minor axis (WGS84)
b_earth = 6356752.3142 # m

# Earth's reverse flattening factor (WGS84)
f_earth = 1/298.257223563 # unitless

# Earth's eccentricity
e_earth = 0.081819221456 # unitless

# Earth's equatorial radius
R_earth = 6378137 # m

# J2 correction term
J2 = 1.08262668e-3 # unitless

# Perturbation constant for Earth's gravity
mu_J2_r_earth_sq = mu_earth * J2 * R_earth**2 # m^5/s^2 

# Earth's rotation rate
omega_earth = 7.292115e-5 # rad/s

# Seconds in a day
sec_in_day = 24 * 60 * 60 # s 

# Seconds in a minute
sec_in_min = 60 # s

# Radians per revolution
rad_per_rev = 2 * np.pi # rad/rev