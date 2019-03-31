#!/usr/bin/env python3

import requests
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from dateutil.parser import parse
from datetime import datetime,timedelta

# GLOBAL VARIABLES
auth_url = 'https://observe.lco.global/api/api-token-auth/'
requ_url = 'https://observe.lco.global/api/userrequests/'
dfmt = '%Y-%m-%d %H:%M:%S'
max_airmass = 2.0
min_lunar_distance = 30.0

###### CHANGE THIS TO YOUR PROPOSAL, USERNAME, PASSWORD!!!!! ########
proposal = 'NOAO2019A-020-TC'
username = 'dummy_user'
password = 'dummy_pass'
#####################################################################

def make_user_request(requests, group_id, typ, ipp, proposal):
    operator = ''
    if len(requests)>1:
        operator = 'MANY'
    elif len(requests)==1:
        operator = 'SINGLE'
    else:
        return(1)

    user_request = {
        'requests': requests,
        'group_id': str(group_id),
        'observation_type': str(typ),
        'operator': operator,
        'ipp_value': float(ipp),
        'proposal': str(proposal)
    }
    return(user_request)

def make_request(location, constraints, target, molecules, windows):
    request = {
        'location': location,
        'constraints': constraints,
        'target': target,
        'molecules': molecules,
        'windows': windows
    }
    return(request)

def make_location(telescope):
    location = {
        'telescope_class': telescope
    }
    return(location)

def make_constraints(max_airmass, min_lunar_distance):
    constraints = {
        'max_airmass': max_airmass,
        'min_lunar_distance': min_lunar_distance
    }
    return(constraints)


# This sets up a target based on name, ra, dec
def make_target(name, ra, dec):
    # parse ra/dec.  LCOGT requires degree formatting
    if ':' in str(ra):
        # Assume hms/dms
        coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    else:
        # Assume decimal degrees
        coord = SkyCoord(ra, dec, unit=(u.deg, u.deg))
    target = {
        'type': 'SIDEREAL',
        'name': name,
        'ra': coord.ra.degree,
        'dec': coord.dec.degree,
        'proper_motion_ra': 0.0,
        'proper_motion_dec': 0.0,
        'parallax': 0.0,
        'coordinate_system': 'ICRS',
        'equinox': 'J2000',
        'epoch': 2000.0,
        'rot_mode': 'VFLOAT',
        'rot_angle': 0.0
    }
    return(target)

# This sets up the instrument parameters for an obs
# Imaging specific
def make_molecule_image(filt, exptime):
    molecule = {
        'type': 'EXPOSE',
        'instrument_name': '1M0-SCICAM-SINISTRO',
        'filter': filt,
        'exposure_time': exptime,
        'exposure_count': 1,
        'bin_x': 1,
        'bin_y': 1,
        'defocus': 0.0,
        'ag_mode': 'OPTIONAL'
    }
    return(molecule)

# Create a window object for LCOGT.  This is defined by a
# start time (earliest the obs can be executed) and a window
# (the amount of time in which obs can be executed)
# Assume duration in days
def make_window(start, duration):
    # Times are formatted like YYYY-MM-DD HH:MM:SS
    startfmt = start.strftime(dfmt)
    endfmt = (start + timedelta(days=duration)).strftime(dfmt)
    window = {
        'start': startfmt,
        'end': endfmt
    }
    return(window)

def post_user_request(user_request, username, password):
    # Authentication
    data = {
        'username': username,
        'password': password
    }
    auth = requests.post(auth_url, data=data).json()
    token = auth['token']

    response = requests.post(requ_url, json=user_request,
        headers={'Authorization': 'Token '+token}).json()
    return(response)

# TEST TEST TEST
# Test, submit 120 second exposures in gr bands to the 1m0 for 2019clr

# TEST TARGET PARAMETERS
obj = {
    'name': '2019clr',
    'ra': 196.4296,
    'dec': 37.6267,
}

# OBSERVING PARAMETERS
obs = {
    'telescope': '1m0',
    'start': datetime.now(),
    'duration': 1,            # in days
    'max_airmass': max_airmass,
    'min_lunar_distance': min_lunar_distance,
    'expose': [{'filt': 'gp', 'exptime': 120}, {'filt': 'rp', 'exptime': 120}]
}

# PACKAGE THESE INTO ONE OBJECT
req = {'obj': obj, 'obs': obs}

# Make objects
target = make_target(req['obj']['name'], req['obj']['ra'], req['obj']['dec'])
window = make_window(req['obs']['start'], req['obs']['duration'])
molecs = []
for exp in req['obs']['expose']:
    molecs.append(make_molecule_image(exp['filt'], exp['exptime']))
locati = make_location(req['obs']['telescope'])
constr = make_constraints(req['obs']['max_airmass'],
    req['obs']['min_lunar_distance'])
request = make_request(locati, constr, target, molecs, [window])
user_request = make_user_request([request], req['obj']['name'],
    'NORMAL', 1.00, proposal)

# Submit request
response = post_user_request(user_request, username, password)

print(response)
