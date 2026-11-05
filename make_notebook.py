"""Build Colonization_History.ipynb, a reduced end-to-end reproduction."""
import json
from pathlib import Path

cells = [
    ("markdown", """# Separating cell-carried from surface-carried colonization history (v1.0)

This notebook reruns a reduced version of every analysis in the package: unit tests,
the calibrated factorial comparison, the mechanism-recovery study and the empirical
re-analysis of Lee *et al.* (2018) from the shipped observables. Reduced ensemble
sizes give noisier numbers than the full runs reported in the manuscript; the full
commands are listed in `README.md`.

Run it from the package directory, or in Colab after uploading and unzipping the
package (first cell)."""),
    ("code", """import os, sys, subprocess, zipfile, glob
if not os.path.exists('model.py'):
    z = sorted(glob.glob('/content/*Colonization_History_Model_v1*.zip'))
    if z:
        zipfile.ZipFile(z[-1]).extractall('/content')
        os.chdir(glob.glob('/content/colonization_history_v1*')[0])
print(os.getcwd())
subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', '-r', 'requirements.txt'])"""),
    ("code", """!python -m unittest test_model -v 2>&1 | tail -25"""),
    ("markdown", """## Calibrated factorial comparison

Uses the shipped calibration (`results/logs/calibration.json`). Writes to `nb_results/`."""),
    ("code", """import shutil, pathlib
out = pathlib.Path('nb_results'); (out/'logs').mkdir(parents=True, exist_ok=True)
shutil.copy('results/logs/calibration.json', out/'logs'/'calibration.json')
!python run_experiments.py --n 100 --check-n 60 --output nb_results | tail -5"""),
    ("code", """from IPython.display import Image, display
for f in ('fig2_matched_baseline_transfer.png', 'fig3_gap_dependence.png'):
    display(Image(filename=f'nb_results/figures/{f}', width=800))"""),
    ("markdown", """## Mechanism recovery (reduced: 400 experiments per hypothesis)

Accuracy is lower than in the full study because matching and training sets are
about seven times smaller."""),
    ("code", """!python recovery.py --n 400 --workers 2 --output nb_results 2>/dev/null | grep -E 'n8\\\\|gaps0-4-16|baseline'"""),
    ("markdown", """## Empirical re-analysis from shipped observables

Ten bootstrap resamples here; the manuscript uses 60."""),
    ("code", """!python empirical/analyze.py 10 nb_results/empirical > /dev/null 2>&1
import json
E = json.load(open('nb_results/empirical/empirical_fits.json'))
for k, v in E['families'].items():
    print(f"{k:28s} k={v['k']:2d} logL={v['loglik']:.1f} AIC={v['aic']:.1f}")
print('adhesiveness shift:', {g: round(x, 3) for g, x in E['families']['selection']['adhesiveness_shift'].items()})
print('LRT selection vs null:', E['lrt_selection_vs_null'])
print('LRT history beyond selection:', E['lrt_history_beyond_selection'])"""),
    ("markdown", """## Full-record and time-matched fits (point estimates only)

Uses the continuous-tracking observables; bootstrap intervals are in the shipped results."""),
    ("code", """import sys; sys.path.insert(0, 'empirical')
import analyze_v2 as V
for mode in ('time', 'occupancy'):
    r = V.time_matched(B=0, mode=mode)
    print(mode, 'offset - head start =', round(r['dll_offset_minus_head'], 1),
          'FC2 offset =', round(r['fits']['offset']['offset']['FC2'], 2))"""),
]

nb = {"cells": [], "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}
for i, (kind, src) in enumerate(cells):
    c = {"cell_type": kind, "id": f"c{i}", "metadata": {}, "source": src.splitlines(keepends=True)}
    if kind == "code":
        c.update({"execution_count": None, "outputs": []})
    nb["cells"].append(c)
Path("Colonization_History.ipynb").write_text(json.dumps(nb, indent=1))
print("notebook written")
