"""Run synthetic calving and melt-mask regressions with a compiled UFEMISM."""
import argparse
import os
from pathlib import Path
import re
import signal
import subprocess

import numpy as np
from netCDF4 import Dataset


def write_input(path, variable, values, times=None):
    axis = np.linspace(-1e6, 1e6, 41)
    with Dataset(path, 'w') as nc:
        for name in ('x', 'y'):
            nc.createDimension(name, len(axis))
            coord = nc.createVariable(name, 'f8', (name,))
            coord[:] = axis
            coord.units = 'm'
        dims = ('y', 'x')
        if times is not None:
            nc.createDimension('time', None)
            time = nc.createVariable('time', 'f8', ('time',))
            time[:] = times
            time.units = 'years'
            dims = ('time',) + dims
        nc.createVariable(variable, 'f8', dims, fill_value=False)[:] = values


def run_case(executable, root, name, kind, fixed, ranks):
    directory = root / name
    directory.mkdir()
    base = Path(__file__).resolve().parents[1] / 'integrated_test_MISMIP_mod_small/config.cfg'
    config = base.read_text()
    settings = {
        'fixed_output_dir_config': repr(str(directory / 'out')),
        'end_time_of_run_config': '1.0',
        'dt_output_config': '0.05', 'dt_output_grid_config': '1.0',
        'dt_climate_config': '0.05', 'dt_ice_max_config': '0.05',
        'dt_ice_min_config': '0.01', 'allow_mesh_updates_config': '.false.',
        'transects_ANT_config': "''",
        'choice_BMB_model_ANT_config': "'prescribed'",
        'choice_BMB_prescribed_ANT_config': "'BMB_no_time'",
        'filename_BMB_prescribed_ANT_config': repr(str(root / 'bmb.nc')),
        'do_asynchronous_BMB_config': '.true.', 'dt_BMB_config': '10.0',
    }
    if fixed:
        settings['choice_ice_integration_method_config'] = "'none'"
    if kind:
        settings.update({
            'do_use_ISMIP_future_shelf_collapse_forcing_config': '.true.',
            'ISMIP_future_shelf_collapse_forcing_filename_config': repr(str(root / 'mask.nc')),
            'shelf_collapse_type_config': repr(kind),
            'retreat_mask_without_time_config': '.false.',
        })
    fields = ['Hi', 'fraction_gr', 'BMB', 'BMB_shelf', 'mask_floating_ice']
    if kind:
        fields.append('retreat_mask')
    settings.update({f'choice_output_field_{i:02d}_config': repr(fields[i-1] if i <= len(fields) else 'none')
                     for i in range(1, 51)})
    extra = []
    for key, value in settings.items():
        config, count = re.subn(r'(?m)^\s*' + key + r'\s*=.*$', f'    {key} = {value}', config)
        if not count:
            extra.append(f'    {key} = {value}')
        elif count != 1:
            raise ValueError(f'Duplicate configuration key: {key}')
    config, count = re.subn(r'(?m)^\s*/\s*$', '\n'.join(extra) + '\n/', config)
    if count != 1:
        raise ValueError('Expected one namelist terminator')
    cfg = directory / 'config.cfg'
    cfg.write_text(config.rstrip() + '\n')
    with (directory / 'run.log').open('w') as log:
        process = subprocess.Popen(['mpirun', '-n', str(ranks), str(executable), str(cfg)],
                                   cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            code = process.wait(timeout=900)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            raise
    if code:
        raise RuntimeError(f'{name} failed with exit {code}; see {directory / "run.log"}')
    with Dataset(directory / 'out/main_output_ANT_00001.nc') as nc:
        result = {field: np.asarray(nc[field][:]) for field in fields + ['V', 'time']}
    if len(result['time']) != 20 or not np.isclose(result['time'][-1], .95):
        raise AssertionError(f'{name}: incomplete output')
    if not np.isfinite(result['Hi']).all():
        raise AssertionError(f'{name}: non-finite thickness')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path, help='New directory for inputs, logs and outputs')
    parser.add_argument('--ranks', type=int, default=2)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    north_grid = np.broadcast_to((np.linspace(-1e6, 1e6, 41) > 0)[:, None], (41, 41))
    times = [0., .2, .3, .6, .7, 1.]
    values = np.array([0., 0., 1., 1., 0., 0.])
    write_input(root / 'mask.nc', 'mask', values[:, None, None] * north_grid, times)
    write_input(root / 'bmb.nc', 'BMB', np.full((41, 41), -5.))
    runs = {name: run_case(executable, root, name, kind, fixed, args.ranks)
            for name, kind, fixed in [('control', None, False), ('calving', 'calving', False),
                                      ('melt', 'BMB', False), ('fixed', 'BMB', True)]}
    control, fixed = runs['control'], runs['fixed']
    for run in runs.values():
        np.testing.assert_array_equal(run['V'], control['V'])
        np.testing.assert_allclose(run['time'], control['time'])
    floating = (control['mask_floating_ice'][0] > 0) & (control['fraction_gr'][0] == 0)
    north = floating & (control['V'][1] > 1e5)
    south = floating & (control['V'][1] < -1e5)
    if not north.any() or not south.any():
        raise AssertionError('Test geometry lacks northern or southern floating ice')
    active = fixed['retreat_mask'][:, north] > .01
    if not active.any() or active.all():
        raise AssertionError('Mask must both activate and deactivate')
    np.testing.assert_allclose(fixed['Hi'], np.broadcast_to(fixed['Hi'][0], fixed['Hi'].shape))
    np.testing.assert_allclose(fixed['BMB'][:, north][active], -100.)
    np.testing.assert_allclose(fixed['BMB'][:, north][~active], -5.)
    np.testing.assert_allclose(fixed['BMB'][:, south], -5.)
    np.testing.assert_allclose(fixed['BMB_shelf'][:, north], -5.)
    expected = np.interp(fixed['time'], times, values)[:, None]
    np.testing.assert_allclose(fixed['retreat_mask'][:, north],
                               np.broadcast_to(expected, fixed['retreat_mask'][:, north].shape), atol=1e-12)
    middle = np.argmin(abs(control['time'] - .5))
    np.testing.assert_allclose(runs['calving']['Hi'][middle, north], 0., atol=1e-8)
    if not np.mean(runs['melt']['Hi'][-1, north]) < np.mean(control['Hi'][-1, north]):
        raise AssertionError('Retreat melt did not thin the selected shelf')
    print('PASS: interpolation, calving, melt cap, deactivation and unselected shelf')


if __name__ == '__main__':
    main()
