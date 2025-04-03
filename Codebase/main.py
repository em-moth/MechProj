import requests
import numpy as np
import json

from utils.utils import get_TLE_data, decode_TLE
from utils.path_utils import obtain_keplerian_elements
from utils.graphing_utils import twoD_kepler_graph, threeD_kepler_graph

data = get_TLE_data(None)

print(data)

TLE = decode_TLE(data['tle'])

k = obtain_keplerian_elements(TLE)

threeD_kepler_graph(k, False)



