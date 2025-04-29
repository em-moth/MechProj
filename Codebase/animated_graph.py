from utils.utils import *
from utils.sgp4 import *
import numpy as np
import matplotlib.pyplot as plt
from utils.graphing_utils import graph_orbit, get_position_array_kepler

from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
data = get_TLE_data(44941)

TLE = decode_TLE(data['tle'])

period = 5
samples = 1000
x = np.empty(samples, dtype=float)
y = np.empty(samples, dtype=float)

p_x = np.empty(samples, dtype=float)
p_y = np.empty(samples, dtype=float)
p_z = np.empty(samples, dtype=float)
t = np.empty(samples, dtype=int)

TEME_pos_i, TEME_vel_i = sgp4(TLE, 0)
gcrs_i = TEME_to_GCRS(TEME_pos_i, TEME_vel_i, TLE['epoch_year'], TLE['epoch_day'], 0)
elements = keplerian_elements(gcrs_i)

#graph_orbit(ax, elements)

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

k_x, k_y, k_z = get_position_array_kepler(elements, t)

tail_length = 10
# Color normalization
norm = Normalize(vmin=0, vmax=t[-1])

# Global scatter holders
scat1 = None
scat2 = None

# Create colorbars manually
sm1 = ScalarMappable(cmap='plasma', norm=norm)
sm1.set_array([])

sm2 = ScalarMappable(cmap='viridis', norm=norm)
sm2.set_array([])


#cbar1 = plt.colorbar(sm1, ax=ax, location='bottom', shrink=0.2, label='Time (min) - SGP4')

#cbar2 = plt.colorbar(sm2, ax=ax, location='bottom', shrink=0.2, label='Time (min) - Kepler orbit')




# Initialize plot
def init():
    ax.set_xlim([-8e3, 8e3])
    ax.set_ylim([-8e3, 8e3])
    ax.set_zlim([-8e3, 8e3])
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title('SGP4 prediction versus Kepler obit of STARLINK-1094')
    return []

# Update plot frame-by-frame
def update(i):
    global scat1, scat2
    if scat1 is not None:
        scat1.remove()
    if scat2 is not None:
        scat2.remove()
    
    start = max(0, i - tail_length)

    # Scatter for Object 1
    scat1 = ax.scatter(
        p_x[start:i+1], p_y[start:i+1], p_z[start:i+1],
        c=t[start:i+1], cmap='plasma', norm=norm, marker='o', s=20
    )
    
    # Scatter for Object 2
    scat2 = ax.scatter(
        k_x[start:i+1], k_y[start:i+1], k_z[start:i+1],
        c=t[start:i+1], cmap='viridis', norm=norm, marker='^', s=20
    )
    
    return [scat1, scat2]

# Create the animation
ani = animation.FuncAnimation(
    fig,
    update,
    frames=samples,
    init_func=init,
    interval=50,  # ms between frames
    blit=False
)

# To save the animation using Pillow as a gif
writer = animation.PillowWriter(fps=15, metadata=dict(artist='Me'),bitrate=1800)
ani.save('gifs\\scatter.gif', writer=writer)


plt.show()