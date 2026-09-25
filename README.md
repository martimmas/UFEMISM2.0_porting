[![License](https://img.shields.io/github/license/UPSY-group/UPSY-models)](LICENSE)
[![Paper](https://img.shields.io/badge/Paper-Published-blue.svg)](https://doi.org/10.5194/gmd-18-3635-2025)
[![Paper](https://img.shields.io/badge/Paper-Published-blue.svg)](https://doi.org/10.5194/egusphere-2026-930)
![GitHub commit activity](https://img.shields.io/github/commit-activity/w/UPSY-group/UPSY-models)
[![Latest Commit](https://img.shields.io/github/last-commit/UPSY-group/UPSY-models)](https://github.com/UPSY-group/UPSY-models/commits)
[![Workflow Status](https://github.com/UPSY-group/UPSY-models/actions/workflows/UFE_test_suite.yml/badge.svg)](https://github.com/UPSY-group/UPSY-models/actions)
[![Workflow Status](https://github.com/UPSY-group/UPSY-models/actions/workflows/UPSY_test_suite.yml/badge.svg)](https://github.com/UPSY-group/UPSY-models/actions)

![UPSY_logo](UPSY_logo.png)

Welcome to the Utrecht Polar SYstem (UPSY) models repo! Here you will find the UPSY modelling toolkit, the
Utrecht Finite Volume Ice-Sheet Model (UFEMISM), and the One-Layer Antarctic Model for Dynamical
Downscaling of Ice–Ocean Exchanges (LADDIE).

See https://github.com/UPSY-group/UPSY-models/wiki/Getting-started for how to set up your model.

### Python tools
Some tools are available to plot model output on its native mesh. To use these, you need
to install [Miniforge3](https://conda-forge.org/download/) and run:
```
conda env create -f environment.yml
conda install -n base -c conda-forge conda-libmamba-solver
conda config --set solver libmamba
conda activate upsy
python -m pip install -e . --no-deps --no-build-isolation
```
To the package can then be loaded by `import upsy`,
You can also try out these commands in the terminal:
```
upsy-diagnose-run rundir
upsy-plot-2dfigure rundir
upsy-plot-3dfigure rundir
```
For additional help, try `upsy-plot-2dfigure -h`

### Prescribed ice-shelf retreat

The retreat-mask forcing, adapted from Franco's `iQ2300_R-LIS` branch, is disabled
by default. Enable it with:

```fortran
do_use_ISMIP_future_shelf_collapse_forcing_config = .true.
ISMIP_future_shelf_collapse_forcing_filename_config = 'retreat.nc'
shelf_collapse_type_config = 'calving' ! or 'BMB'
retreat_mask_without_time_config = .false.
```

Supply a NetCDF `mask` in [0,1], without missing values, on an x/y grid, lon/lat
grid or model mesh. Transient masks need strictly increasing times in model
years. Frames are interpolated linearly; endpoints are held outside the time
range. Values above 0.01 select retreat, sampled at the climate timestep.

`calving` removes selected floating ice. `BMB` prescribes −400 m/yr shelf melt
before subgrid weighting and the global melt cap (100 m/yr by default); the
background melt is restored when selection ends. Retreat BMB cannot be combined
with `inverted` or `prescribed_fixed` BMB. Request `retreat_mask` as an output
field only when the forcing is enabled.

For static masks, `retreat_mask_applied_only_to_open_ocean_config = .true.`
restricts selection to initially ice-free ocean. Reuse the generated
`retreat_mask_open_ocean_reference_<REGION>.nc` through
`retreat_mask_open_ocean_reference_filename_config` when restarting. Use
`{region}` in the configured filename for multiple regions. Remeshing preserves
this initial reference rather than deriving it from the evolving geometry.

See [retreat-mask tests](automated_testing/UFEMISM/retreat_mask/README.md).
