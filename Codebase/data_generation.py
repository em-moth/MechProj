from utils.utils import *
from utils.sgp4 import *
import numpy as np
import matplotlib.pyplot as plt
from utils.graphing_utils import graph_orbit

from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

ax = plt.figure().add_subplot(projection='3d')

data = get_TLE_data(None)

TLE = decode_TLE(data['tle'])

period = 40
samples = 200
x = np.empty(samples, dtype=float)
y = np.empty(samples, dtype=float)

p_x = np.empty(samples, dtype=float)
p_y = np.empty(samples, dtype=float)
p_z = np.empty(samples, dtype=float)
t = np.empty(samples, dtype=int)

TEME_pos_i, TEME_vel_i = sgp4(TLE, 0)
gcrs_i = TEME_to_GCRS(TEME_pos_i, TEME_vel_i, TLE['epoch_year'], TLE['epoch_day'], 0)
elements = keplerian_elements(gcrs_i)

graph_orbit(ax, elements)

for i in range(0, samples):
    t[i] = i* period
    TEME_pos, TEME_vel = sgp4(TLE, i * period)
    gcrs = TEME_to_GCRS(TEME_pos, TEME_vel, TLE['epoch_year'], TLE['epoch_day'], i * period)
    s = (gcrs.cartesian.xyz.to(u.km).value.flatten())
    #ax.scatter(s[0], s[1], s[2], c='r', marker='o')
    p_x[i] = s[0]
    p_y[i] = s[1]
    p_z[i] = s[2] 
    new_elements = propagate_elements(elements, i * period)
    k = kepler_elements_to_coords(new_elements)

    dist = (np.sqrt((s[0] - k[0]) ** 2 + (s[1] - k[1]) ** 2 + (s[2] - k[2]) ** 2))

    x[i] = i * period
    y[i] = dist

colors = t
sc = ax.scatter3D(p_x,p_y,p_z, c=colors, cmap='viridis', marker='o', label='SGP4 prediction')
plt.colorbar(sc, ax=ax, label='Time (min)', location='bottom', shrink=0.35)

ax.set(xlabel='X position (km)', ylabel='Y position (km)', zlabel='Z position (km)')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))

# Add the legend
#ax.legend(loc=3)



plt.show()
# Create the scatter plot
plt.scatter(x, y)

# Add labels and title
plt.xlabel("Time (min)")
plt.ylabel("Distance (km)")
plt.title("Distance between SGP4 propogated orbit and Kepler propogated orbit of the ISS over time")



# Show the plot
plt.show()