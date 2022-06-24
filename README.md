# Atmospheric moisture tracking

This repository contains the source code of the WAM2layers moisture tracking
code originally developed by Ruud van der Ent. It can be used to determine where
precipitation originally evaporated (backtracking), or where evaporated moisture
eventually ends up (forward tracking).

## How to use

The model workflow consists of three steps:

1. Preprocess the input data. This step converts raw input data to a format that
   is understood by WAM-2layers (described in detail below). Preprocessing
   scripts are readily available for several commonly used datasets. See
   `preprocess_era5.py` as an example. Some useful utility functions are
   available in the file `preprocessing.py`.
2. Run the tracking code. The script `backtrack.py` performs backtracking for a
   selected region. Forward tracking is currently missing from this codebase,
   but we are planning to reinstate it soon(ish).
3. Optionally postprocess the output data. Some legacy code is available here:
   - `Con_E_Recyc_Output_monthly_mat.py`: aggregate daily data as monthly means
     and store in convenient format.
   - `Hor_Fluxes_Output_mat.py`: convert intermediate output (moisture fluxes over the grid cells) from step 1
     preprocessing to monthly data.


## Raw data format (ERA-5)
For ERA-5 the following variables are read in and converted to the flux and moisture variables used by the backtracking script:

    - u: zonal wind in m/s
    - v: meridonal wind in m/s
    - q: specific humidity in kg kg-1
    - sp: surface pressure in Pa
    - evap: evaporation in m (accumulated hourly)
    - cp: convective precipitation in m (accumulated hourly)
    - lsp: large scale precipitation in m (accumulated hourly)
If you want to apply WAM2layers to a different data set, please use `preprocess_era5.py` as a basis and change accordingly to your data set.

## Backtracking input data format
The backtracking code makes several assumptions about the incoming data. We use
the term pre-processing for anything that's done to convert raw input data to
the format that is supported by WAM-2layers. The final input data into the
tracking should have the following characteristics.

- Data should be stored in netcdf files, one file per day.
- The file should contain the following variables:
   - fa_e_upper: eastward moisture transport in the upper layer
   - fa_n_upper: northward moisture transport in the upper layer
   - fa_e_lower: eastward moisture transport in the lower layer
   - fa_n_lower: northward moisture transport in the lower layer
   - w_upper: total water in the upper layer
   - w_lower: total water in the lower layer
   - fa_vert: vertical transport of water between the layers
   - evap: evaporation from the surface
   - precip: precipiation
- All variables should be in units of m3
- Flux variables are in between the state variables, e.g. for hourly data the
  state variables (w_upper and w_lower) come at 00, 01, ... 00, and the fluxes
  come at 00:30, 01:30, .. 23:30. For the state variables, the midnight edges
  for one day and the next are included in both files, so there is some
  duplication there.
- Precipitation and evaporation should both be positive.

Here is an example of a preprocessed netCDF file. Note that the latitude,
longitude, and time may vary for your data.

```
Dimensions:     (time: 95, lat: 121, lon: 321, time2: 96)
Dimensions without coordinates: time, lat, lon, time2
Data variables:
    fa_e_upper  (time, lat, lon) float64 ...
    fa_n_upper  (time, lat, lon) float64 ...
    fa_e_lower  (time, lat, lon) float64 ...
    fa_n_lower  (time, lat, lon) float64 ...
    w_upper     (time2, lat, lon) float64 ...
    w_lower     (time2, lat, lon) float64 ...
    fa_vert     (time, lat, lon) float64 ...
    evap        (time, lat, lon) float32 ...
    precip      (time, lat, lon) float32 ...
```

## Other versions

This is the official codebase for the WAM-2layers moisture tracking model as of
18/03/2022, but there are still several other versions around:

- [Original Python code by Ruud van der Ent](https://github.com/ruudvdent/WAM2layersPython)
- [Adapted version by Imme Benedict](https://github.com/Imme1992/moisture_tracking_mississippi)

## Reuse and acknowledgement
To be completed.

We are actively developing the code at the moment, so it may be subject to
change. We encourage anyone who is interested in re-using the code to get in
touch. We may be able to help.

If you use the code for a publication, please cite it using its DOI (TODO)
and/or the following paper: [Contrasting roles of interception and transpiration
in the hydrological cycle - Part 2: Moisture
recycling](https://doi.org/10.5194/esd-5-471-2014)
