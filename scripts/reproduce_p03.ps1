# Run from repository root. A clean checkout has no results/ directory.
# Windows Python 3.10; no merge/deployment actions.
$ErrorActionPreference = 'Stop'
function Invoke-CheckedPython {
    param([string[]] $PythonArgs)
    & .venv/Scripts/python.exe @PythonArgs
    if ($LASTEXITCODE -ne 0) { throw "Python command failed: $PythonArgs" }
}
Invoke-CheckedPython -PythonArgs @('scripts/fetch_p03_data.py')
Invoke-CheckedPython -PythonArgs @('scripts/validate_sampler.py')
Invoke-CheckedPython -PythonArgs @('scripts/reproduce_juliet.py')
Invoke-CheckedPython -PythonArgs @('scripts/check_tutorial_likelihood.py')
Invoke-CheckedPython -PythonArgs @('scripts/validate_tutorial_direct.py')
Invoke-CheckedPython -PythonArgs @('scripts/validate_tutorial_direct.py', '--models', '2', '--nlive', '2000', '--sample', 'unif', '--seed-base', '19141', '--tag', 'global')
Invoke-CheckedPython -PythonArgs @('scripts/run_p03.py', '--models', 'A', 'B', 'C', 'N', '--seeds', '303', '73344')
Invoke-CheckedPython -PythonArgs @('scripts/run_p03.py', '--variant', 'broad', '--models', 'A', 'B', 'C', '--seeds', '303', '73344')
Invoke-CheckedPython -PythonArgs @('scripts/run_p03.py', '--variant', 'amplitude', '--models', 'A', 'B', 'C', '--seeds', '303')
Invoke-CheckedPython -PythonArgs @('scripts/run_p03.py', '--models', 'C', '--nlive', '1000', '--seeds', '1303', '2303')
Invoke-CheckedPython -PythonArgs @('scripts/plot_p03.py')
Invoke-CheckedPython -PythonArgs @('scripts/summarize_p03.py')
