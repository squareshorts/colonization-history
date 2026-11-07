$ErrorActionPreference = 'Stop'
python -m pip install -r requirements.txt
if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed' }
python -m unittest discover -s tests -v
if ($LASTEXITCODE -ne 0) { throw 'Unit tests failed' }
python src/run_capability_sweep.py --reps 48
if ($LASTEXITCODE -ne 0) { throw 'Capability sweep failed' }
python -c "import sys;sys.path.insert(0,'src');from run_robustness import run;run(10)"
if ($LASTEXITCODE -ne 0) { throw 'Robustness analysis failed' }
python src/oscillatory_signatures.py
if ($LASTEXITCODE -ne 0) { throw 'Oscillator analysis failed' }
python src/run_surface_conditioning.py
if ($LASTEXITCODE -ne 0) { throw 'Surface-conditioning control failed' }
python src/make_results.py
if ($LASTEXITCODE -ne 0) { throw 'Figure/summary generation failed' }
Write-Output 'Reproduction complete; see REPORT.md, results, and figures.'
