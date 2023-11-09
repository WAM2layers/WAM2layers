### About how test data is prepared

The [test data](./era5/) is already included in this directory, which can be used to run the [tests](../test_workflow.py).

It is prepared through the following steps:
- Download the required ERA5 fields from CDS <br>
For model levels the following variables are needed:
  - u, v, q on selected levels
  - tp, e, sp, and tcw at the surface
(see the example downloading script for more details about CDS query of ERA5 fields in `scripts/download_era5_ml.py`. Note that in the tests we use hourly data of corresponding ERA5 fields of the date 2022-08-31.)

- Slice data to have a minimum sample set for testing and generate `source_region.nc`.

### License
Test data included here is generated using Copernicus Climate Change Service information and for more information about licensing, please check the [Licence Agreement](https://cds.climate.copernicus.eu/cdsapp/#!/terms/licence-to-use-copernicus-products) for Copernicus Products.
