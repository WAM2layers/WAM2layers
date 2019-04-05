# moisture_tracking_mississippi
Change of WAM-2layers model to work with climate model data on five pressure levels

The scripts in this folder are an adaption of the WAM-2layers model developed by Ruud van der Ent from Delf University (https://github.com/ruudvdent/WAM2layersPython).
An explanation on the code can be found in his dissertation (Van der Ent, R. J. (2014), A new view on the hydrological cycle over continents, Ph.D. thesis, 96 pp, Delft University of Technology, Delft. 
http://dx.doi.org/10.4233/uuid:0ab824ee-6956-4cc3-b530-3245ab4f32be.).

The code was originally written to work with re-analysis model data from ERA-Interim.
I have adapted the code such that it can run with 'regular output' of global climate models, namely atmospheric output on pressure levels.

I have used data from the climate model EC-Earth (Hazeleger et al., 2010, 2014), which was available at 5 pressure levels in the atmosphere,
namely 850, 700, 500, 300 and 200 hPa.
We use (logarithmic) surface pressure to eliminate levels which are located below the surface.
We apply a spline interpolation on the moisture fluxes (u*q and v*q) to correct for information that we miss close to the surface (such as low level jets).
We apply a linear interpolation on the specific humidity vertical profiles.
The atmospheric data has a time-step of 6 hours, the surface data a time-step of 3 hours.
