from utils.utils import *
from utils.sgp4 import *
import numpy as np
import matplotlib.pyplot as plt

data = get_TLE_data(None)

TLE = decode_TLE(data['tle'])

period = 5
samples = 60
x = np.empty(samples, dtype=float)
y = np.empty(samples, dtype=float)

TEME_pos_i, TEME_vel_i = sgp4(TLE, 0)
gcrs_i = TEME_to_GCRS(TEME_pos_i, TEME_vel_i, TLE['epoch_year'], TLE['epoch_day'], 0)
elements = keplerian_elements(gcrs_i)

for i in range(0, samples):
    TEME_pos, TEME_vel = sgp4(TLE, i * period)
    gcrs = TEME_to_GCRS(TEME_pos, TEME_vel, TLE['epoch_year'], TLE['epoch_day'], i * period)
    s = (gcrs.cartesian.xyz.to(u.km).value.flatten())
    
    new_elements = propagate_elements(elements, i * period)
    k = kepler_elements_to_coords(new_elements)

    dist = (np.sqrt((s[0] - k[0]) ** 2 + (s[1] - k[1]) ** 2 + (s[2] - k[2]) ** 2))

    x[i] = i * period
    y[i] = dist


# Create the scatter plot
plt.scatter(x, y)

# Add labels and title
plt.xlabel("Time (min)")
plt.ylabel("Distance (km)")
plt.title("Distance between SGP4 propogated orbit and Kepler propogated orbit of the ISS over time")



# Show the plot
plt.show()