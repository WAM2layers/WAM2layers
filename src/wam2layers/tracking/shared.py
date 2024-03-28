import pandas as pd
import xarray as xr


def initialize_time(config, direction="forward"):
    dt = pd.Timedelta(seconds=config.timestep)

    if direction == "forward":
        t0 = pd.Timestamp(config.tracking_start_date)
        th = t0 + dt / 2  # th is "half" time, i.e. between t0 and t1
        t1 = t0 + dt
    elif direction == "backward":
        t1 = pd.Timestamp(config.tracking_end_date)
        th = t1 - dt / 2  # th is "half" time, i.e. between t0 and t1
        t0 = t1 - dt
    else:
        raise ValueError("Direction should be forward or backward")

    return t0, th, t1, dt


def initialize_tagging_region(bbox, lat, lon):
    """Build a new xr.DataArray with ones inside and zeros outside bbox."""
    lat_in_bbox = (bbox.south <= lat) & (lat <= bbox.north)
    unrolled_lon_in_bbox = (bbox.west <= lon) & (lon <= bbox.east)
    rolled_lon_in_bbox = (bbox.west <= lon) | (lon <= bbox.east)
    rolled_domain = bbox.west > bbox.east
    lon_in_bbox = rolled_lon_in_bbox if rolled_domain else unrolled_lon_in_bbox
    tagging_region = xr.where((lat_in_bbox & lon_in_bbox), 1, 0).values
    return tagging_region
