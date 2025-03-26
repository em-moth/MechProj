import requests
import numpy as np
import json

from Codebase.utils.utils import get_TLE_data, decode_TLE

data = get_TLE_data(None)

print(data)

print(decode_TLE(data['tle']))



