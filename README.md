# Atmospheric moisture tracking

This repository contains the source code of the WAM2layers moisture tracking
code originally developed by Ruud van der Ent. It can be used to determine where
precipitation originally evaporated (backtracking), or where evaporated moisture
eventually ends up (forward tracking).


## Solved equations
Considering the schematic illustration below, forward- and backward tracking can be formulated respectively as:

$$
\begin{align}
s_{t+1}=s_t &+F_{we,w}\ast s_{t,w}-F_{we,e}\ast s_t & s_{t-1}=s_t &+F_{we,e}\ast s_{t,e}-F_{we,w}\ast s_t \\
            &+F_{ew,e}\ast s_{t,e}-F_{ew,w}\ast s_t &             &+F_{ew,w}\ast s_{t,w}-F_{ew,e}\ast s_t \\
            &+F_{ns,n}\ast s_{t,n}-F_{ns,s}\ast s_t &             &+F_{ns,s}\ast s_{t,s}-F_{ns,n}\ast s_t \\
            &+F_{sn,s}\ast s_{t,s}-F_{sn,n}\ast s_t &             &+F_{sn,n}\ast s_{t,n}-F_{sn,s}\ast s_t \\
            &+F_{ul,u}\ast s_{t,u}-F_{ul,l}\ast s_t &             &+F_{ul,l}\ast s_{t,l}-F_{ul,u}\ast s_t \\
            &+F_{lu,l}\ast s_{t,l}-F_{lu,u}\ast s_t &             &+F_{lu,u}\ast s_{t,u}-F_{lu,l}\ast s_t \\
            &+E-P                                   &             &-E+P
\end{align}
$$

where $F$ are (total) moisture fluxes and $s$ represents the amount of tracked relative to total moisture in the grid cells. 
The up/down fluxes are not illustrated but they follow the same systematic as the north/south and east/west fluxes. 
Note that all fluxes are by definition positive; this is needed because moisture flux is scaled with the relative amount of tracked moisture 
in the grid cells where it originates from.

In WAM2Layers code, these equations are solved for one upper and one lower layer (such that only two of the 4 vertical transport terms are relevant for each layer). The evaporation term is used for the lower layer only, while the precipitation contribution is distributed across the two layers.

![image](https://user-images.githubusercontent.com/17080502/183941265-116cfb6d-3b03-4602-a750-9abe7c7cc554.png)

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
    - evap: evaporation in m (accumulated hourly): Be aware evaporation in ERA5 has a negative sign
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
