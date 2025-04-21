from .constants import API_KEY_NY2GO
import requests
import json
from astropy.coordinates import ITRS, GCRS, CartesianRepresentation, CartesianDifferential
from astropy.time import Time
from astropy import units as u
import numpy as np
import spiceypy as spice
import os
import skyfield.sgp4lib as sgp4lib
import erfa

# Get absolute path to data folder
DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))


# Load kernels
spice.furnsh(os.path.join(DATA_DIR, 'latest_leapseconds.tls'))
spice.furnsh(os.path.join(DATA_DIR, 'gm_de431.tpc'))


ny2go_api_key = API_KEY_NY2GO

# Function to get data from the API
def get_data_dict(request_url):
    response = requests.get(f"https://api.n2yo.com/rest/v1/satellite{request_url}&apiKey={ny2go_api_key}")

    if response.status_code == 200:
        try:
            data_dict = response.json()
            # Process the dictionary
            return data_dict
        except ValueError:
            raise Exception("Response is not in JSON format")
    else:
        raise Exception(f"Request failed with status code {response.status_code}")

# Function to get the TLE data for a satellite
def get_TLE_data(norad_id):
    if (norad_id == None):
        return {'info': {'satid': 25544, 'satname': 'SPACE STATION', 'transactionscount': 3}, 'tle': '1 25544U 98067A   25085.46276703  .00029905  00000-0  52819-3 0  9998\r\n2 25544  51.6371 356.7528 0003651  53.2985 306.8339 15.50092829502334'}
    return get_data_dict(f"/tle/{norad_id}")

#Function to decode a TLE
def decode_TLE(tle):
    data_dict = {}
    lines = tle.split("\r\n")
    if len(lines) != 2:
        raise Exception("Invalid TLE format")
    
    # Extract the two lines
    line_1 = lines[0]
    line_2 = lines[1]

    # Extract the data from the first line
    data_dict["satellite_number"] = int(line_1[2:7])
    data_dict["classification"] = line_1[7]

    data_dict["international_designator"] = line_1[9:17]
    data_dict["epoch_year"] = int(line_1[18:20])
    data_dict["epoch_day"] = float(line_1[20:32])
    
    data_dict["first_derivative_of_mean_motion"]=float(line_1[33:43])
    data_dict["second_derivative_of_mean_motion"]=sci_to_float(line_1[44:52])
    
    data_dict["drag_term"]=sci_to_float(line_1[53:61])

    data_dict["ephemeris_type"] = int(line_1[62])

    data_dict["element_set_number"] = int(line_1[64:68])

    data_dict["line_1_checksum"] = int(line_1[68])

    #Extract data from the second line

    data_dict["inclination"] = float(line_2[8:16])
    data_dict["right_ascension_of_ascending_node"] = float(line_2[17:25])
    data_dict["eccentricity"] = float(f"0.{line_2[26:33]}")
    data_dict["argument_of_perigee"] = float(line_2[34:42])

    data_dict["mean_anomaly"] = float(line_2[43:51])
    data_dict["mean_motion"] = float(line_2[52:63])

    data_dict["rev_number_at_epoch"] = int(line_2[63:68])

    data_dict["line_2_checksum"] = int(line_2[68])

    return data_dict

#convert scientific notation string to float
def sci_to_float(s):
    significand_str = s[:-2]
    exponent_str = s[-2:]
    significand = float(significand_str) / 1e5  # Add this division!
    exponent = int(exponent_str)
    return significand * (10 ** exponent)

def TEME_to_GCRS(r_teme, v_teme, epoch_year, epoch_day, delta_t_mins):
    """
    Convert position and velocity vectors from TEME to GCRS coordinates.

    Parameters
    ----------
    r_teme : numpy.ndarray
        TEME position vector (3 elements) in kilometers.
    v_teme : numpy.ndarray
        TEME velocity vector (3 elements) in kilometers per second.
    epoch_year : int
        Two-digit year (YY format, >=57 implies 19YY, <57 implies 20YY).
    epoch_day : float
        Day of year with fractional days.
    delta_t_mins : float
        Time correction in minutes to add to epoch_day.

    Returns
    -------
    astropy.coordinates.GCRS
        GCRS coordinate object containing position and velocity.
    """
    # Time conversion
    full_year = 1900 + epoch_year if epoch_year >= 57 else 2000 + epoch_year
    epoch_day += delta_t_mins / 1440  # Convert minutes to fractional days

    # Create Astropy Time object
    time_str = (
        f"{full_year}:{int(epoch_day):03d}:"
        f"{int((epoch_day % 1) * 24):02d}:"
        f"{int(((epoch_day % 1) * 24 % 1) * 60):02d}:"
        f"{((epoch_day % 1) * 86400) % 60:09.6f}"
    )
    t = Time(time_str, format="yday", scale="utc").tt
    # Get precession-nutation matrix using IAU 2006/2000A model
    # Convert time to ERFA format (JD1, JD2)
    jd1, jd2 = t.jd1, t.jd2
    rnpb = erfa.pnm06a(jd1, jd2)

    # Rotate vectors from TEME to GCRS
    position_gcrs = rnpb @ r_teme  # Matrix multiplication
    velocity_gcrs = rnpb @ v_teme

    # Create GCRS object with proper units
    return GCRS(
        CartesianRepresentation(
            x=position_gcrs[0] * u.km,
            y=position_gcrs[1] * u.km,
            z=position_gcrs[2] * u.km,
            differentials=CartesianDifferential(
                d_x=velocity_gcrs[0] * u.km/u.s,
                d_y=velocity_gcrs[1] * u.km/u.s,
                d_z=velocity_gcrs[2] * u.km/u.s
            )
        ),
        obstime=Time('J2000', scale='tt')
    )

    

def keplerian_elements(gcrs):
    # Ensure velocity exists
    if not hasattr(gcrs, 'velocity') or gcrs.velocity is None:
        raise ValueError("GCRS object must include velocity information")

    # Extract position (km) and velocity (km/s) in GCRS frame (≈J2000)
    r_km = gcrs.cartesian.xyz.to(u.km).value.flatten()
    v_kms = gcrs.velocity.d_xyz.to(u.km/u.s).value.flatten()
    state = np.concatenate([r_km, v_kms]).astype(np.float64)
    
    # Get Earth's GM using NAIF ID 399 (Earth)
    # Requires prior kernel load: spice.furnsh("de421.bsp")
    mu = spice.bodvrd('399', 'GM', 1)[1][0]  # km³/s²
    
    # Convert GCRS obstime (TT scale) to UTC for SPICE
    # SPICE requires UTC for str2et (assumes UTC input)
    t_utc = gcrs.obstime.utc
    et = spice.str2et(t_utc.isot)  # Correct epoch in TDB seconds
    
    # Compute osculating elements (state in J2000 frame)
    elements = spice.oscelt(state, et, mu)
    
    # Calculate semi-major axis (km)
    rp = elements[0]  # Periapsis radius (km)
    ecc = elements[1]
    a_km = rp / (1 - ecc)  # Valid only for eccentricity < 1
    
    return (
        a_km,         # Semi-major axis [km]
        ecc,          # Eccentricity
        elements[2],  # Inclination [rad]
        elements[3],  # RAAN [rad]
        elements[4],  # Arg of periapsis [rad]
        elements[5],  # Mean anomaly at epoch [rad]
        et            # Epoch (seconds past J2000 TDB)
    )

def true_anomaly_from_elements(elements, tolerance=1e-12, max_iter=1000):
    """
    Calculate true anomaly from Keplerian elements tuple.
    
    Parameters:
    elements : tuple (a_km, e, i, Ω, ω, M0, et_initial)
        - M0: Mean anomaly at epoch [rad]
    tolerance : float, optional
        Convergence tolerance for solving Kepler's equation
    max_iter : int, optional
        Maximum number of iterations for Newton-Raphson
        
    Returns:
    float: True anomaly in radians [0, 2π)
    """
    # Extract needed parameters
    a_km, e, _, _, _, M0, _ = elements
    
    # Handle circular case
    if e < 1e-12:
        return M0 % (2 * np.pi)
    
    # Solve Kepler's equation E - e*sin(E) = M using Newton-Raphson
    E = M0  # Initial guess
    for _ in range(max_iter):
        delta = (E - e * np.sin(E) - M0) / (1 - e * np.cos(E))
        E -= delta
        if np.abs(delta) < tolerance:
            break
    
    # Calculate true anomaly with quadrant awareness
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E/2),
                          np.sqrt(1 - e) * np.cos(E/2))
    
    return nu % (2 * np.pi)
    


def propagate_elements(elements, delta_t_mins):
    """
    Propagate Keplerian elements forward in time using two-body dynamics.
    
    Parameters:
    elements : tuple (a_km, e, i, Ω, ω, M0, et_initial)
        - a_km: Semi-major axis [km]
        - e: Eccentricity
        - i: Inclination [rad]
        - Ω: RAAN [rad]
        - ω: Argument of perigee [rad]
        - M0: Mean anomaly at epoch [rad]
        - et_initial: Initial epoch [seconds since J2000]
    delta_t_mins : float
        Propagation time in minutes
        
    Returns:
    tuple: Updated elements in same format (a_km, e, i, Ω, ω, M_new, et_new)
    """
    # Unpack elements
    a_km, e, i, Ω, ω, M0, et_initial = elements
    
    # Earth's gravitational parameter (km³/s²)
    mu = 398600.4418
    
    # Convert time and calculate mean motion
    delta_t_sec = delta_t_mins * 60
    n = np.sqrt(mu / a_km**3)  # Mean motion [rad/s]
    
    # Update mean anomaly and epoch
    M_new = (M0 + n * delta_t_sec) % (2 * np.pi)
    et_new = et_initial + delta_t_sec
    
    return (a_km, e, i, Ω, ω, M_new, et_new)

#note that all arguments should be in radians
def perifocal_to_equatorial(x, y, z, inclination, arg_periapsis, LAAN):
    #Omega = LAN
    #omega = arg_periapsis
    #i = inclination
    rotation_matrix = np.matrix([
        [
            np.cos(LAAN) * np.cos(arg_periapsis) - np.sin(LAAN) * np.cos(inclination) * np.sin(arg_periapsis),
            -1 *np.cos(LAAN) * np.sin(arg_periapsis) - np.sin(LAAN) * np.cos(inclination) * np.cos(arg_periapsis),
            np.sin(LAAN) * np.sin(inclination)
        ],
        [
            np.sin(LAAN) * np.cos(arg_periapsis) + np.cos(LAAN) * np.cos(inclination) * np.sin(arg_periapsis),
            -1 * np.sin(LAAN) * np.sin(arg_periapsis) + np.cos(LAAN) * np.cos(inclination) * np.cos(arg_periapsis),
            -1* np.cos(LAAN)*np.sin(inclination)
        ],
        [
            np.sin(inclination) * np.sin(arg_periapsis),
            np.sin(inclination) * np.cos(arg_periapsis),
            np.cos(inclination)
        ]
    ])

    equatorial_eq = [0, 0, 0]
    equatorial_eq[0] = rotation_matrix[0, 0] * x + rotation_matrix[0, 1] * y + rotation_matrix[0, 2] * z
    equatorial_eq[1] = rotation_matrix[1, 0] * x + rotation_matrix[1, 1] * y + rotation_matrix[1, 2] * z
    equatorial_eq[2] = rotation_matrix[2, 0] * x + rotation_matrix[2, 1] * y + rotation_matrix[2, 2] * z
    return equatorial_eq

def kepler_elements_to_coords(kepler_elements):
    a, e, i, raan, argp, nu, _ = kepler_elements

    v = true_anomaly_from_elements(kepler_elements)

    #Perifocal coordinate system https://en.wikipedia.org/wiki/Perifocal_coordinate_system
    r = (a * (1 - e**2)) / (1 + e * np.cos(v))

    x = r * np.cos(v)
    y = r * np.sin(v)
    z = 0

    equatorial_eq_v = perifocal_to_equatorial(x, y, z, i, argp, raan)
    x_rotated = equatorial_eq_v[0]
    y_rotated = equatorial_eq_v[1]
    z_rotated = equatorial_eq_v[2]

    return (x_rotated, y_rotated, z_rotated)