### About how test data is prepared

The [test data](./era5/) is already included in this directory, which can be used to run the [tests](../test_workflow.py).

It is prepared through the following steps:
- Download the required ERA5 fields from CDS <br>
For model levels the following variables are needed:
  - u, v, q on selected levels
  - tp, e, sp, and tcw at the surface
(see the example downloading script for more details about CDS qurey of ERA5 fields in `scripts/download_era5_ml.py`)

- Slice data to have a minimum sample set for testing and generate `source_region.nc`.
