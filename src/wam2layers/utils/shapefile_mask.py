"""Utils to create a tagging region mask from a shapefile on the fly."""
from collections import namedtuple
from pathlib import Path
from typing import Union

import numpy as np
import shapefile
import shapely.geometry
import xarray as xr

Resolution = namedtuple("Resolution", ["lat", "lon"])


def ensure_monotonic(
    obj: Union[xr.DataArray, xr.Dataset],
) -> Union[xr.DataArray, xr.Dataset]:
    """Ensure that both latitude and longitude are both monotonically increasing."""
    for coord in ("latitude", "longitude"):
        if not obj.indexes[coord].is_monotonic_increasing:
            obj = obj.sortby(coord)
        if not obj.indexes[coord].is_unique:
            obj = obj.drop_duplicates(coord)
    return obj


def load_shapefile(file: Union[str, Path]):
    """ "Load the first polygon from a shapefile."""
    sf = shapefile.Reader(file)
    return shapely.geometry.shape(sf.shape(0))


def infer_resolution(ds: xr.Dataset) -> Resolution:
    return Resolution(
        lat=abs(float(ds["latitude"].diff(dim="latitude", n=1).median())),
        lon=abs(float(ds["longitude"].diff(dim="longitude", n=1).median())),
    )


def generate_boxes(
    da_loc: xr.DataArray,
    resolution: Resolution,
) -> list[shapely.Polygon]:
    """Generate boxes for each gridcell."""
    boxes = []
    for loc in da_loc:
        x = loc["longitude"].to_numpy()
        y = loc["latitude"].to_numpy()
        box = shapely.box(
            x - 0.5 * resolution.lon,
            y - 0.5 * resolution.lat,
            x + 0.5 * resolution.lon,
            y + 0.5 * resolution.lat,
        )
        boxes.append(box)
    return boxes


def create_mask(ds: xr.Dataset, shape: Union[str, Path]):
    """Create a tagging region mask from a shapefile."""
    poly = load_shapefile(shape)
    minx, miny, maxx, maxy = poly.bounds
    resolution = infer_resolution(ds)

    subset = ensure_monotonic(ds).sel(
        latitude=slice(miny - resolution.lat, maxy + resolution.lat),
        longitude=slice(minx - resolution.lon, maxx + resolution.lon),
    )
    stack = subset.stack(loc=("longitude", "latitude"))
    boxes = generate_boxes(stack["loc"], resolution)

    weights = []
    for box in boxes:
        weights.append(shapely.intersection(box, poly).area / box.area)

    subset_mask = xr.DataArray(
        np.array(weights).reshape((subset["longitude"].size, subset["latitude"].size)),
        {"longitude": subset["longitude"], "latitude": subset["latitude"]},
    )
    subset_mask.name = "tagging_region"
    mask = subset_mask.interp(
        {"latitude": ds["latitude"].data, "longitude": ds["longitude"].data},
        method="nearest",
    )
    dim_order = [dim for dim in ds.dims if dim in ("latitude", "longitude")]
    mask = mask.transpose(*dim_order)  # make sure the dims are correctly ordered.
    return mask.fillna(0)
