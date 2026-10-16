# PHANTOM-PLANET P01

Research feasibility prototype. Synthetic stellar radial-velocity signals, candidate-period analysis and an optional juliet fitting script. This is **not** a validated planet detector.

## Run

```bash
pip install -r requirements.txt
python run_p01.py --seeds 10
python -m unittest discover -s tests -v
```

See the P01 review branch for implementation details. Synthetic outputs are reproducible and are not astronomical discoveries.
