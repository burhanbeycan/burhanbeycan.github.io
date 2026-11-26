"""Check the monotone response surface against measured conditions and untouched holdouts."""
import importlib.util
import json
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parents[1]
web = root / 'dist' if (root / 'dist').is_dir() else root
spec = importlib.util.spec_from_file_location('hpam_shape', root / 'scripts/train-constrained-hpam.py')
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)
x, y = model.load(web / 'viscosity-inputs/hpam-akbari.csv')
actual = model.benchmark(x, y)
saved = json.loads((web / 'family-inputs/hpam-shape.json').read_text())
assert len(saved['rows']) == len(actual['rows']) == 15
assert len({r['group'] for r in actual['rows']}) == 9
for left,right in zip(saved['rows'],actual['rows']):
    assert left['id'] == right['id']
    np.testing.assert_allclose([left['candidate'],left['baseline']], [right['candidate'],right['baseline']], rtol=1e-6)
for point in actual['rows']:
    error = abs(point['candidate'] - point['observed'])
    assert np.isfinite(error)
z = np.log(x / model.REFERENCE)
assert np.min(model.derivative_constraints(z) @ np.asarray(actual['coef'])) >= -1e-7
hold = (x[:,0] == 2000) & (x[:,1] == 12.5)
adjusted = y.copy()
adjusted[hold] *= 3
before = [r['candidate'] for r in actual['rows'] if r['group']=='2000 ppm · 12.5 g/L NaCl']
after = [r['candidate'] for r in model.benchmark(x,adjusted)['rows'] if r['group']=='2000 ppm · 12.5 g/L NaCl']
assert len(before) == 3
np.testing.assert_allclose(before,after)
