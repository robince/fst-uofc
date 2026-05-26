# fst-uofc

Python 3.13 bindings for the General Fourier Family Transform implementation in
this repository.

## Development

```bash
uv sync
.venv/bin/python -m pytest -q
```

The Python package imports as `pygft`.

```python
import numpy as np
import pygft

signal = np.sin(2 * np.pi * 8 * np.linspace(0, 1, 256))
transform = pygft.gft1d(signal)
```

Interactive Python examples are in `examples/` and use `# %%` cell markers for
VS Code/Jupyter-style execution.
