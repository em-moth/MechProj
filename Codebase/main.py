import requests
import numpy as np
import json

from utils.utils import get_TLE_data, decode_TLE, TEME_to_GCRS, get_true_anomaly, keplerian_elements
from utils.graphing_utils import twoD_kepler_graph, threeD_kepler_graph
from utils.sgp4 import sgp4
from astropy import units as u

data = get_TLE_data(None)
#data_css = get_TLE_data(54216)

print(data)

TLE = decode_TLE(data['tle'])
print(TLE)

TEME_pos, TEME_vel = sgp4(TLE, 0)
gcrs = TEME_to_GCRS(TEME_pos, TEME_vel, TLE['epoch_year'], TLE['epoch_day'], 0)
print(gcrs.cartesian.xyz.to(u.km).value.flatten())
print(gcrs.velocity.d_xyz.to(u.km/u.s).value.flatten())
print(gcrs)
elements = keplerian_elements(gcrs)
threeD_kepler_graph([elements], True)



