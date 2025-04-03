from .constants import API_KEY_NY2GO
import requests
import json

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
def sci_to_float(str):
    return float(str[:-2]) * pow(10,int(str[-2:]))