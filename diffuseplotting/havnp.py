import numpy as np
from numba import autojit

# we only want longitude in order to configure the plothealpix_map module.
@autojit
def haversine(lon1, lat1, lon2, lat2, onEarth=False):
    # Calculate the great circle distance between two points
    # on the earth (specified in decimal degrees)

    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    if onEarth:
        r = 6371000  # Radius of earth in meters.
        return c * r # Return physical distance between lat/lon points on Earth
    else:
        return c* 180.0/np.pi # Return angular distance in degrees.
