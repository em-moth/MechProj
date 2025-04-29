from utils.utils import decode_TLE
from astropy.time import Time

def read_TLE_data(filename, start, end):
    data = []
    try:
        with open(filename, 'r') as file:
            while True:
                line1 = file.readline()
                line2 = file.readline()

                if not line1 and not line2:
                    break  # End of file

                line1 = line1.strip()
                line2 = line2.strip()

                tle = decode_TLE(f"{line1}\r\n{line2}")
                epoch_year = tle['epoch_year']
                epoch_day = tle['epoch_day']
                full_year = 1900 + epoch_year if epoch_year >= 57 else 2000 + epoch_year
                # Create Astropy Time object
                time_str = (
                    f"{full_year}:{int(epoch_day):03d}:"
                    f"{int((epoch_day % 1) * 24):02d}:"
                    f"{int(((epoch_day % 1) * 24 % 1) * 60):02d}:"
                    f"{((epoch_day % 1) * 86400) % 60:09.6f}"
                )
                t = Time(time_str, format="yday", scale="utc")
                if start.utc < t < end.utc:
                    tle['time'] = t
                    data.append(tle)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    return data


