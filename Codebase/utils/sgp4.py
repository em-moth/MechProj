import numpy as np



def sgp4(ELE, t):
    deg2rad = np.pi / 180
    #extract elements from TLE
    n_0 = ELE["mean_motion"] * 2 * np.pi / 1440.0
    e_0 = ELE["eccentricity"]
    i_0 = ELE["inclination"] * deg2rad 
    M_0 = ELE["mean_anomaly"] * deg2rad
    omega_0 = ELE["argument_of_perigee"] * deg2rad
    Omega_0 = ELE["right_ascension_of_ascending_node"] * deg2rad
    nd_0 = ELE["first_derivative_of_mean_motion"]
    ndd_0 = ELE["second_derivative_of_mean_motion"]
    Bstar = ELE["drag_term"]
    
    #obtain some constants

    #radius of the earth
    ER = 6378.137

    #gravitational constant
    G = 6.67430 * (10 ** -11)

    #mass of the earth
    M = 5.972 * (10 ** 24)

    a_E_ER = 1.0  # 1 Earth radius = 6378.135 km
    XKMPER = 6378.135  # km/ER (conversion factor)

    #k_e = np.sqrt(G * M)
    k_e = 0.0743669161
    #zonal harmonics
    J_2 = 1.08262668 * (10 ** -3)
    J_3 = -2.5324105 * (10 ** -6)
    J_4 = -1.6198976 * (10 ** -6)

    #more constants
    k_2 = (1/2) * J_2 * (a_E_ER ** 2)

    k_4 = (-3/8) * J_4 * (a_E_ER ** 4)

    A_30 = -1 * J_3 * (a_E_ER ** 3)

    q_0 = 120.0 / XKMPER  # Critical fix: ~0.01882 ER


    s = 1.01222928 * a_E_ER

    #Calculations
    a_1 = (k_e / n_0) ** (2/3)

    delta_1 = (3/2) * (k_2 / (a_1 ** 2)) * (3 * (np.cos(i_0) ** 2) - 1) / ((1 - e_0 ** 2) ** (3/2))

    a_0 = a_1 * (1 - (delta_1/3) - (delta_1 **2) - (134/81) * (delta_1 ** 3))

    delta_0 = (3/2) * (k_2 / (a_0 ** 2)) * ((3 * (np.cos(i_0) ** 2) - 1)/((1 - e_0 ** 2) ** (3/2)))

    naa_0 = n_0 / (1 + delta_0)

    aaa_0 = a_0 / (1- delta_0)

    r_perigee_km = (aaa_0*(1.0 - e_0) - a_E_ER) * XKMPER  # km
    s = a_E_ER
    QOS= (120.0 - 78.0)/XKMPER**4  # Default Q0MS2T
    
    if r_perigee_km < 98.0:
        s = 20.0/XKMPER + a_E_ER
        QOS = ((120.0 - 20.0)/XKMPER)**4
    elif r_perigee_km < 156.0:
        s = (r_perigee_km - 78.0)/XKMPER + a_E_ER
        QOS = ((120.0 - (r_perigee_km - 78.0))/XKMPER)**4


    theta = np.cos(i_0)

    xi = 1/(aaa_0 - s)

    beta_0 = (1 - e_0 **2) ** (1/2)

    eta = aaa_0 * e_0 * xi

    

    C2 = QOS * (xi ** 4) * naa_0 * ((1-eta**2) ** (-7/2)) * ( aaa_0*(1 + (3* eta ** 2)/2 + 4 * e_0 * eta + e_0 * (eta **3) + (3/2) * ((k_2 * xi)/(1-eta ** 2)) * (-1/2 + (3/2) * theta ** 2) * (8 + 24 * eta **2 + 3 * eta ** 4)) )

    C1 = Bstar * C2

    C3 = QOS * xi ** 5 * A_30 * naa_0 * a_E_ER * np.sin(i_0) / (k_2 * e_0)

    C4 = 2 * naa_0 * QOS * xi ** 4 * aaa_0 * beta_0 ** 2 * ((1 - eta **2) ** (-7/2)) *( (2 * eta * (1 + e_0 * eta) + e_0/2 + (eta ** 3)/2) - ((2 * k_2 * xi)/(aaa_0 * (1-eta ** 2))) * (3 * (1 - 3 * theta ** 2) * (1 + (3/2) * eta ** 2 - 2 * e_0 * eta - (1/2) * e_0 * eta - (1/2) * e_0 * eta ** 3) + (3/4) * (1-theta ** 2)*(2 * eta **2 - e_0 * eta - e_0 * eta ** 3)*np.cos(2*omega_0)) )

    C5 = 2 * QOS * xi ** 4 * aaa_0 * beta_0 ** 2 * (1- eta **2) ** (-7/2) * (1 + (11/4)*eta * (eta + e_0)+e_0 * eta**3)

    D2 = 4 * aaa_0 * xi* C1 ** 2

    D3 = (4/3) * aaa_0 * xi ** 2 * (17 * aaa_0 + s) * C1 ** 3

    D4 = (2/3) * aaa_0 * xi ** 3 *(221 * aaa_0 + 31 * s) * C1 **4

    M_DF = M_0 + (1 + (3 * k_2 * (-1 + 3 * theta ** 2))/(2 * aaa_0 ** 2 * beta_0 ** 4) + (3 * k_2 ** 2 * (7 - 114 * theta ** 2 + 395 * theta ** 4))/(16* aaa_0 **4 * beta_0 ** 8) + (5*k_4 *(3 - 36*theta ** 2 + 49 * theta **4))/(4*aaa_0 ** 4 * beta_0 ** 8)) * naa_0 * t

    omega_DF = omega_0 + ( -1 * (3 * k_2 * (1 - 5 * theta ** 2))/(2*aaa_0 ** 2 * beta_0 ** 4) + (3 * k_2 **2 * (7-114 * theta **2 + 395 * theta **4))/(16 * aaa_0 ** 4 * beta_0 ** 8) + (5*k_4 * (3 - 36 * theta **2 + 49 * theta ** 4))/(4 * aaa_0 ** 4 * beta_0 ** 8)) * naa_0 * t

    Omega_DF = Omega_0 + ( -1 * (3*k_2 * theta)/(aaa_0 ** 2 * beta_0 ** 4) + (3 * k_2 **2 *(4*theta-19*theta**3))/(2*aaa_0**4*beta_0**8) + (5*k_4*theta * (3-7*theta **2))/(2 * aaa_0**4 * beta_0 ** 8) ) * naa_0 * t

    
    if r_perigee_km < 220:
        M_P = M_DF
        omega = omega_DF
        Omega = Omega_DF - (21/2) * ((naa_0 * k_2 * theta)/(aaa_0 ** 2 * beta_0 ** 2)) * C1 * t ** 2
        e = e_0 - Bstar * C4 * t

        a = aaa_0 * (1 - C1 * t) ** 2
        L = M_P + omega + Omega + naa_0 * (3/2) * C1 * t ** 2
    else:
        deltaomega = Bstar * C3 * np.cos(omega_0) * t

        deltaM = (-2/3) * QOS * Bstar * xi ** 4 * (a_E_ER/(e_0 * eta)) * ((1 + eta * np.cos(M_DF)) ** 3 - (1 + eta * np.cos(M_0)) ** 3)

        M_P = M_DF + deltaomega + deltaM
        
        omega = omega_DF - deltaomega - deltaM

        Omega = Omega_DF - (21/2) * ((naa_0 * k_2 * theta)/(aaa_0 ** 2 * beta_0 ** 2)) * C1 * t ** 2

        e = e_0 - Bstar * C4 * t - Bstar * C5 * (np.sin(M_P) - np.sin(M_0))

        a = aaa_0 * (1 - C1 * t - D2 * t ** 2 - D3 * t **3 - D4 * t **4) ** 2

        L = M_P + omega + Omega + naa_0 * ( (3/2) * C1 * t ** 2 + (D2 + 2 * C1 **2) * t ** 3 + (1/4) * (3 * D3 + 12 * C1 *D2 + 10 * C1  ** 3) * t **4 + (1/5) * (3 * D4 + 12 * C1 * D3 + 6*D2 **2 + 30 * C1 **2 * D2 + 15 * C1 ** 4) * t ** 5 )
    
    beta = np.sqrt(1 - e ** 2)

    n = k_e/(a ** (3/2))

    a_xN = e * np.cos(omega)

    L_L = (A_30 * np.sin(i_0))/(8 * k_2 * a * beta ** 2) * e * np.cos(omega) * ((3 + 5 * theta)/ 1 + theta)

    a_yNL = (A_30 * np.sin(i_0))/(4 * k_2 * a * beta  ** 2)

    L_T = L + L_L

    a_yN = e * np.sin(omega) + a_yNL

    U = L_T - Omega

    EW = U
    EW_next = EW + (U - a_yN * np.cos(EW) + a_xN * np.sin(EW) - EW)/(-1 * a_yN * np.sin(EW) - a_xN * np.cos(EW) + 1)
    
    while np.abs(EW_next - EW) > 1e-12:
        EW = EW_next
        EW_next = EW + (U - a_yN * np.cos(EW) + a_xN * np.sin(EW) - EW)/(-1 * a_yN * np.sin(EW) - a_xN * np.cos(EW) + 1)
    
    EW = EW_next

    ecosE = a_xN * np.cos(EW) + a_yN * np.sin(EW)

    esinE = a_xN * np.sin(EW) - a_yN * np.cos(EW)

    e_L = (a_xN ** 2 + a_yN ** 2) ** (1/2)

    p_L = a * (1-e_L ** 2)

    r = a * (1 - ecosE)

    rdot = k_e * (np.sqrt(a)/r) * esinE

    rfdot = k_e * (np.sqrt(p_L))/r

    cosu = (a/r) * (np.cos(EW) - a_xN + (a_yN * esinE)/(1 + np.sqrt(1 - e_L ** 2)))

    sinu = (a/r) * (np.sin(EW) - a_yN - (a_xN * esinE)/(1 + np.sqrt(1-e_L ** 2)))

    u = np.arctan2(sinu, cosu)

    Deltar = (k_2)/(2*p_L) * (1 - theta **2) * np.cos(2 * u)

    Deltau = -1 * (k_2)/(4 * p_L **2) * (7 * theta ** 2 -1) * np.sin(2 * u)
    DeltaOmega = (3 * k_2 * theta)/(2 * p_L ** 2) * np.sin(2 * u)
    Deltai = (3 * k_2 * theta)/(2 * p_L ** 2) * np.sin(2 * u)

    Deltardot =-1 * (k_2 * n)/(p_L) * (1 - theta ** 2) * np.sin(2 * u)

    Deltarfdot = (k_2 * n)/(p_L)* ((1-theta ** 2) * np.cos(2 * u)-(3/2) * (1 - 3 *theta ** 2))

    r_k = r * (1 - (3/2) * k_2 * (np.sqrt(1-e_L**2))/(p_L ** 2) * (3 * theta ** 2 - 1)) + Deltar

    u_k = u + Deltau

    Omega_k = Omega + DeltaOmega

    i_k = i_0 + Deltai

    rdot_k = rdot + Deltardot

    rfdot_k = rfdot + Deltarfdot

    M_vec =np.array([[-1 * np.sin(Omega_k) * np.cos(i_k)],
    [np.cos(Omega_k) * np.cos(i_k)],
    [np.sin(i_k)]])

    N_vec = np.array([[np.cos(Omega_k)],
    [np.sin(Omega_k)],
    [0]])

    U_vec = M_vec * np.sin(u_k) + N_vec * np.cos(u_k)
    V_vec = M_vec * np.cos(u_k) - N_vec * np.sin(u_k)

    position = r_k * U_vec * XKMPER
    velocity = (rdot_k * U_vec + rfdot_k * V_vec) * XKMPER / 60

    return (position, velocity)



def test_sgp4():
    # Example TLE-like orbital elements (values are approximated for testing)
    #ELE = {'satellite_number': 25544, 'classification': 'U', 'international_designator': '98067A  ', 'epoch_year': 25, 'epoch_day': 85.46276703, 'first_derivative_of_mean_motion': 0.00029905, 'second_derivative_of_mean_motion': 0.0, 'drag_term': 0.0005281900000000001, 'ephemeris_type': 0, 'element_set_number': 999, 'line_1_checksum': 8, 'inclination': 51.6371, 'right_ascension_of_ascending_node': 356.7528, 'eccentricity': 0.0003651, 'argument_of_perigee': 53.2985, 'mean_anomaly': 306.8339, 'mean_motion': 15.50092829, 'rev_number_at_epoch': 50233, 'line_2_checksum': 4}
    ELE = {
        "mean_motion": 16.05824518,          # Revolutions per day (line 2, field 7)
        "eccentricity": 0.0086731,           # Line 2, field 4 (0.0086731)
        "inclination": 72.8435,              # Degrees (line 2, field 2)
        "mean_anomaly": 110.5714,            # Degrees (line 2, field 6)
        "argument_of_perigee": 52.6988,      # Degrees (line 2, field 5)
        "right_ascension_of_ascending_node": 115.9689,  # Degrees (line 2, field 3)
        "first_derivative_of_mean_motion": 0.00146188,  # 2 × 0.00073094 (line 1, field 4)
        "second_derivative_of_mean_motion": 0.00083064, # 6 × 0.00013844 (line 1, field 5 "13844-3")
        "drag_term": 6.6816e-05              # 0.66816e-4 from "66816-4" (line 1, field 6)
    }
    # Time since epoch in minutes
    for t_minutes in range(0, 1440+360, 360):

        # Run SGP4
        position_vector, velocity_vector = sgp4(ELE, t_minutes)

        # Output result
        print(f"{t_minutes} - Position vector (km):")
        print(position_vector.flatten())
        print(f"{t_minutes} - Velocity vector (km/min):")
        print(velocity_vector.flatten())

    

# Run the test
test_sgp4()


