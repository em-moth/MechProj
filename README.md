# MechProj
Codebase for Cooper Union Physics 1: Mechanics final project

Command line notes

To activate venv: .\venv\Scripts\activate

To deactivate venv: deactivate

Example TLE dict:

{'info': 
    {'satid': 25544, 'satname': 'SPACE STATION', 'transactionscount': 3}, 'tle': '1 25544U 98067A   25085.46276703  .00029905  00000-0  52819-3 0  9998\r\n2 25544  51.6371 356.7528 0003651  53.2985 306.8339 15.50092829502334'}

Example processed TLE dict:
{
    'satellite_number': 25544, 
    'classification': 'U', 
    'international_designator': '98067A  ', 
    'epoch_year': 25, 
    'epoch_day': 85.46276703, 
    'first_derivative_of_mean_motion': 0.00029905, 
    'second_derivative_of_mean_motion': 0.0, 
    'drag_term': 52.819, 
    'ephemeris_type': 0, 
    'element_set_number': 999, 
    'line_1_checksum': 8, 
    'inclination': 51.6371, 
    'right_ascension_of_ascending_node': 356.7528, 
    'eccentricity': 0.0003651, 
    'argument_of_perigee': 53.2985, 
    'mean_anomaly': 306.8339, 
    'mean_motion': 15.50092829, 
    'rev_number_at_epoch': 50233, 
    'line_2_checksum': 4
}