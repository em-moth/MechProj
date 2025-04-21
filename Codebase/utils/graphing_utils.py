import numpy as np
import matplotlib.pyplot as plt
from utils.utils import perifocal_to_equatorial


def twoD_kepler_graph(kepler_elements):
    # Unpack the tuple into named variables
    a, e, i, raan, argp, v, _ = kepler_elements

    E = np.linspace(0, 2 * np.pi, 1000)
    # Compute r using the polar equation of an ellipse
    r = (a * (1 - e**2)) / (1 + e * np.cos(E))

    # Set up the polar plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    # Plot the orbit
    ax.plot(E, r, color="blue")

    v_r = ((a * (1 - e**2)) / (1 + e * np.cos(v)))
    ax.plot(v, v_r, marker='o', markersize=10, color='red')


    earth_radius = 6356

    r_circle = np.full_like(E, earth_radius)  # Circle with constant radius

    # Plot the circular reference orbit
    ax.plot(E, r_circle, linestyle="dashed", color="gray")
    

    ax.grid(True)
    
    
    ax.legend()
    # Show the plot
    plt.show()


def threeD_kepler_graph(kepler_elements_list, show_earth=False):
    
    earth_radius = 6356

    ax = plt.figure().add_subplot(projection='3d')


    #Show vernal equinox
    for k in kepler_elements_list:
        graph_orbit(ax, k)
    #Show earth
    if show_earth:
        u = np.linspace(0, 2 * np.pi, 10)
        v = np.linspace(0, np.pi, 10)
        x_earth = earth_radius * np.outer(np.cos(u), np.sin(v))
        y_earth = earth_radius * np.outer(np.sin(u), np.sin(v))
        z_earth = earth_radius * np.outer(np.ones(np.size(u)), np.cos(v))

        # Plot the surface
        ax.plot_wireframe(x_earth, y_earth, z_earth, color='green')
        

    plt.axis('equal')
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.legend()

    plt.show()

def graph_orbit(ax, kepler_elements):
    E = np.linspace(0, 2 * np.pi, 1000)
    # Unpack the tuple into named variables
    a, e, i, raan, argp, v, _ = kepler_elements

    #Perifocal coordinate system https://en.wikipedia.org/wiki/Perifocal_coordinate_system
    r = (a * (1 - e**2)) / (1 + e * np.cos(E))
    r_v = (a * (1 - e**2)) / (1 + e * np.cos(v))

    x = r * np.cos(E)
    y = r * np.sin(E)
    z = 0

    x_v = r_v * np.cos(v)
    y_v = r_v * np.sin(v)
    z_v = 0

    equatorial_eq_v = perifocal_to_equatorial(x_v, y_v, z_v, i, argp, raan)
    x_v_rotated = equatorial_eq_v[0]
    y_v_rotated = equatorial_eq_v[1]
    z_v_rotated = equatorial_eq_v[2]

    print(x_v_rotated)
    print(y_v_rotated)
    print(z_v_rotated)

    equatorial_eq = perifocal_to_equatorial(x, y, z, i, argp, raan)

    x_rotated = equatorial_eq[0]
    y_rotated = equatorial_eq[1]
    z_rotated = equatorial_eq[2]

    ax.scatter(x_v_rotated, y_v_rotated, z_v_rotated, c='r', marker='o') 
    ax.plot(x_rotated, y_rotated, z_rotated)
    
    



#[ 4647.8511381  -3306.80826047 -3699.45697685]