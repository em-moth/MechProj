import math


#Get traditional Keplerian elements from TLE set (passed as dict)

def obtain_keplerian_elements(TLE):
    #LAN is the longitude of the ascending node
    LAN = TLE['right_ascension_of_ascending_node'] * (math.pi / 180)

    #Argument of periapsis
    AP = TLE['argument_of_perigee'] * (math.pi / 180)

    #Convert mean motion to radians per second
    MM_rad_s = TLE['mean_motion'] * (1 / 86400) * 2 * math.pi

    #Obtain the period of the orbit
    #period = MM_rad_s * 2 * math.pi 

    #grav_parameter_earth = 3.986004418e14

    #Find the semimajor axis
    #semimajor_axis = (((period ** 2) * grav_parameter_earth)/(4 * (math.pi ** 2))) ** (1/3)
    # Convert mean motion from revs/day to rad/s
    MM_rad_s = TLE['mean_motion'] * (2 * math.pi / 86400)

    # Gravitational parameter for Earth (m^3/s^2)
    grav_parameter_earth = 3.986004418e14

    # Semi-major axis calculation (in meters)
    semimajor_axis= (grav_parameter_earth / (MM_rad_s ** 2)) ** (1 / 3)

    print(semimajor_axis)
    #Obtain true anomaly

    kepler_elements = {}
    kepler_elements["Omega"] = LAN
    kepler_elements["omega"] = AP
    kepler_elements["a"] = semimajor_axis
    kepler_elements["e"] = TLE["eccentricity"]
    kepler_elements["i"] = TLE["inclination"] * (math.pi / 180)
    return kepler_elements


