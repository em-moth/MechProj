from utils.read_data import read_TLE_data
from astropy.time import Time
import numpy as np
from utils.utils import *
from utils.sgp4 import *
import matplotlib.pyplot as plt
from utils.graphing_utils import graph_orbit
from astropy import units as u
import csv
import sys

from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D

if (len(sys.argv) < 4):
    sys.exit()

filename_in = sys.argv[1]
filename_out = sys.argv[2]
sat_name = sys.argv[3]

ax = plt.figure().add_subplot(projection='3d')

start = Time("2000-01-01 00:00:05", scale="utc")
end = Time("2000-04-20 00:00:05", scale="utc")
data = (read_TLE_data(filename_in, start, end))

TLE_i = data[0]
t_i = TLE_i['time']

samples = len(data)

kep_to_actual = np.empty(samples, dtype=float)
sgp_to_actual = np.empty(samples, dtype=float)

p_x = np.empty(samples, dtype=float)
p_y = np.empty(samples, dtype=float)
p_z = np.empty(samples, dtype=float)
t = np.empty(samples, dtype=int)

c_x = np.empty(samples, dtype=float)
c_y = np.empty(samples, dtype=float)
c_z = np.empty(samples, dtype=float)

TEME_pos_i, TEME_vel_i = sgp4(TLE_i, 0)
gcrs_i = TEME_to_GCRS(TEME_pos_i, TEME_vel_i, TLE_i['epoch_year'], TLE_i['epoch_day'], 0)
elements = keplerian_elements(gcrs_i)

graph_orbit(ax, elements)

for i in range(len(data)):
    #get time since start
    TLE_current = data[i]
    t_current = TLE_current['time']
    delta_t = (t_current - t_i).to_value(u.min)

    #calculate position from current TLE
    TEME_pos, TEME_vel = sgp4(TLE_current, 0)
    gcrs = TEME_to_GCRS(TEME_pos, TEME_vel, TLE_current['epoch_year'], TLE_current['epoch_day'], 0)
    sg = (gcrs.cartesian.xyz.to(u.km).value.flatten())
    c_x[i] = sg[0]
    c_y[i] = sg[1]
    c_z[i] = sg[2] 

    #Calculate SGP4 propagation from initial TLE
    t[i] = delta_t
    TEME_pos, TEME_vel = sgp4(TLE_i, delta_t)
    gcrs = TEME_to_GCRS(TEME_pos, TEME_vel, TLE_i['epoch_year'], TLE_i['epoch_day'], delta_t)
    s = (gcrs.cartesian.xyz.to(u.km).value.flatten())
    p_x[i] = s[0]
    p_y[i] = s[1]
    p_z[i] = s[2] 

    #calculate Kepler propogation from initial elements
    new_elements = propagate_elements(elements, delta_t)
    k = kepler_elements_to_coords(new_elements)

    dist_k = (np.sqrt((sg[0] - k[0]) ** 2 + (sg[1] - k[1]) ** 2 + (sg[2] - k[2]) ** 2))
    dist_sgp = (np.sqrt((sg[0] - s[0]) ** 2 + (sg[1] - s[1]) ** 2 + (sg[2] - s[2]) ** 2))
    
    kep_to_actual[i] = dist_k
    sgp_to_actual[i] = dist_sgp

colors = t
sc = ax.scatter3D(p_x,p_y,p_z, c=colors, cmap='viridis', marker='o')
sc_ = ax.scatter3D(c_x,c_y,c_z, c=colors, cmap='viridis', marker='^')
plt.colorbar(sc, ax=ax, label='Time')

plt.show()
# Create the scatter plot
plt.scatter(t, kep_to_actual)
plt.scatter(t, sgp_to_actual)

# Add labels and title
plt.xlabel("Time (min)")
plt.ylabel("Distance (km)")
plt.title(f"Distance between {sat_name} position and SGP4 and Kepler propogated positions over time")
plt.legend(['Kepler prediction', 'SGP4 prediction'])

#do trendlines
z1 = np.polyfit(t, kep_to_actual, 1)
p1 = np.poly1d(z1)
plt.plot(t, p1(t), "r--")

equation = f"y = {z1[0]:.2f}x + {z1[1]:.2f}"
print("Kepler equation:")
print(equation)

z2 = np.polyfit(t, sgp_to_actual, 1)
p2 = np.poly1d(z2)
plt.plot(t, p2(t), "g--")

equation = f"y = {z2[0]:.2f}x + {z2[1]:.2f}"
print("SGP4 equation:")
print(equation)



# Show the plot
plt.show()


with open(filename_out, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile,lineterminator="\n")
    # writing the fields
    csvwriter.writerow(['time', 'Kepler dist', 'SGP4 dist'])
    for i in range(len(t)):
        csvwriter.writerow([t[i], kep_to_actual[i], sgp_to_actual[i]])