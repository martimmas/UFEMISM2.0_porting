# Retreat-mask tests

The Fortran unit tests run through `UFEMISM_unit_test_program` and cover NetCDF
time types, interpolation, endpoint holding, remapping and open-ocean references.

For a short MISMIP_mod integration test, build UFEMISM normally, then run:

```sh
python3 automated_testing/UFEMISM/retreat_mask/test_retreat_mask.py \
    /path/to/UFEMISM_program /path/to/new-test-directory --ranks 2
```

Requires NumPy, netCDF4 and `mpirun`. Four one-year runs check calving, melt
activation/deactivation, the melt cap and preservation of unselected shelf melt.
Inputs are synthetic; configurations, logs and outputs stay in the supplied
directory, which must not already exist. Each run has a 15-minute timeout.
