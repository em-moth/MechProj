import requests
import numpy as np
import json

from utils.utils import get_TLE_data, decode_TLE
from utils.path_utils import obtain_keplerian_elements
from utils.graphing_utils import twoD_kepler_graph, threeD_kepler_graph

data = get_TLE_data(63383)
#data_css = get_TLE_data(54216)

print(data)

TLE = decode_TLE(data['tle'])
#TLE2 = decode_TLE(data_css['tle'])

k = obtain_keplerian_elements(TLE)
#k_css = obtain_keplerian_elements(TLE2)

twoD_kepler_graph(k)
threeD_kepler_graph([k], True)



