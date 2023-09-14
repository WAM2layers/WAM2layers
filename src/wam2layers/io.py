from dataclasses import dataclass
from functools import partial
from typing import Callable, Optional
import xarray as xr
from pandas import Timestamp


def _unpack(t: Timestamp) -> dict[str, int]:
    """Convert pd.Timestamp to dict with individual datatime components."""
    return dict(
        year=t.year, month=t.month, day=t.day, hour=t.hour, minute=t.minute
    )


# TODO cache this function?
def raw_load(filename_template: str, variable: str, time: Optional[Timestamp]=None, **kwargs) -> xr.DataArray:
    """Open a dataset that contains given variable and optionally extract time"""
    if time is not None:
        kwargs.update(**_unpack(time))
    filename = filename_template.format(variable=variable, **kwargs)
    return xr.open_dataset(filename)[variable]


def maybe_interp(loader: Callable, t: Optional[Timestamp]=None, freq: Optional[str]=None) -> xr.DataArray:
    """Wraps a data loader to deal with time.

    Arguments:
        loader: function that can load data for a given timestamp
        t: timestamp
        freq: frequency of the input data used to find adjacent points.

    Returns:
        Data array at time t. If the time is not readily available from input
        data, the data will be interpolated between adjacent points.
    """

    # No time selection
    if t is None:
        return loader(time=t)

    # Timestamp readily available
    t0 = t.floor(freq)  # type: ignore typenarrowing for optional doesn't work
    da0 = loader(time=t0).sel(time=t0, drop=True)
    if t == t0:
        return da0

    # Interpolate data to timestamp
    t1 = t.ceil(freq)  # type: ignore typenarrowing for optional doesn't work
    da1 = loader(time=t1).sel(time=t1, drop=True)
    return da0 + (t - t0) / (t1 - t0) * (da1 - da0)


def load_input(config, variable, time):
    options = {'variable': variable, 'dataset': config.dataset}
    loader = partial(raw_load, filename_template=config.data_reference['input'], **options)
    # TODO rename levels and select subset of levels?
    return maybe_interp(loader, t=time, freq=config.data_frequency['input'])


def load_preprocessed(config, variable, time):
    options = {'variable': variable, 'dataset': config.dataset}
    loader = partial(raw_load, filename_template=config.data_reference["preprocessed"], **options)
    return maybe_interp(loader, t=time, freq=config.data_frequency['preprocessed'])


def load_output(config, variable, time):
    options = {'variable': variable, 'dataset': config.dataset}
    loader = partial(raw_load, filename_template=config.data_reference["output"], **options)
    return maybe_interp(loader, t=time, freq=config.data_frequency['output'])


def load_tracking_region(config, time=None):
    options = {'variable': 'tracking_region', 'dataset': config.dataset}
    loader = partial(raw_load, filename_template=config.data_reference["tracking_region"], **options)
    return maybe_interp(loader, t=time, freq=config.data_frequency['tracking_region'])


def load_source_region(config, time=None):
    options = {'variable': 'source_region', 'dataset': config.dataset}
    loader = partial(raw_load, filename_template=config.data_reference["source_region"], **options)
    return maybe_interp(loader, t=time, freq=config.data_frequency['source_region'])


# TODO: add write functions as well


if __name__ == "__main__":

    @dataclass
    class DataConfig:
        dataset: str
        data_reference: dict[str, str]
        data_frequency: dict[str, Optional[str]]


    config = DataConfig(
        dataset = 'ERA5',
        data_reference = dict(
            input = "~/wam2layers/input_data/{dataset}/{year}/{month:02}/{day:02}/{dataset}_{year}-{month:02d}-{day:02d}_{variable}.nc",
            preprocessed = "~/wam2layers/preprocessed_data/{dataset}/{year}/{month:02}/{day:02}/{dataset}_{variable}.nc",
            output = "~/wam2layers/output_data/{identifier}/{year:04}-{month:02}-{day:02}T{hour:02}-{minute:02}.nc",
            source_region = "~/wam2layers/auxiliary_data/{dataset}/source_region_global.nc",
            tracking_region = "~/wam2layers/auxiliary_data/{dataset}/tracking_region_global.nc",
        ),
        data_frequency = dict(
            input = "1h",
            preprocessed = "1h",
            output = "1h",
            source_region = None,
            tracking_region = None,
        )
    )

    t = Timestamp('20210701T0100')
    th = Timestamp('20210701T0130')

    u_ml = load_input(config, 'u_ml', time=th)
    q_ml = load_input(config, 'q_ml', time=t)

    precip = load_preprocessed(config, 'precip', time=th)
    evap = load_preprocessed(config, 'evap', time=th)

    source_region = load_source_region(config, time=t)
    tracking_region = load_tracking_region(config, time=t)

    load_output(config, 's_track', time=t)
    load_output(config, 'e_track', time=t)
