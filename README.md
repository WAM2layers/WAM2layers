# Atmospheric moisture tracking

This repository contains the source code of the WAM2layers moisture tracking
code originally developed by Ruud van der Ent. It can be used to determine where
precipitation originally evaporated (backtracking), or where evaporated moisture
eventually ends up (forward tracking).

## How to use

The current version of the code is tailored to EC-Earth climate model output.
It assumes an input directory with the following files:

- EVAP_200201_NH.nc (evaporation at the surface)
- LNSP_200201_NH.nc (logarithm of surface pressure)
- Q_200201_NH.nc (specific humidity at pressure levels)
- Q2M_200201_NH.nc (specific humidity at the surface)
- TP_200201_NH.nc (convective precipitation at the surface)
- U_200201_NH.nc (eastward wind field at pressure levels)
- U10_200201_NH.nc (eastward wind at the surface)
- V_200201_NH.nc (northward wind field at pressure levels)
- V10_200201_NH.nc (northward wind at the surface)
- landseamask_ECearth_T799.nc (land-sea mask)
- mask_Miss_ECEarth.npy (mask to designate area of interest)

The workflow is as follows:

1. Preprocess the input data. The script `Fluxes_and_States_mat.py` converts the
   data to 2 layers and calculates the moisture fluxes between these layers and
   between all grid cells. It also does time interpolation to make sure the CFL
   criterion is not violated during the model run. Output is stored in
   intermediate (matlat) files.
2. Run the tracking code. The script `Backtrack_savedaily_mat.py` performs
   backtracking. Forward tracking is currently missing from this codebase.
3. Postprocess the output data.
   - `Con_E_Recyc_Output_monthly_mat.py`: aggregate daily data as monthly means
     and store in convenient format.
   - `Hor_Fluxes_Output_mat.py`: convert intermediate output from step 1
     preprocessing to monthly data.

## Other versions

This is the official codebase for the WAM2Layers moisture tracking model as of
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
