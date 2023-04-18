This configuration can be used for testing wam2layers.

From the root of the repo, run 

$ wam2layers backtrack era5 sample_data/sample_config.yaml 
$ wam2layers backtrack sample_data/sample_config.yaml

Output will be stored in temp dir (works on linux).

Reference output is stored in this same folder, this can be used to
verify that the data hasn't changed.

This data is used in tests/test_cli.py
