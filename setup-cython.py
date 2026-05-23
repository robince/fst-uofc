"""Compatibility entrypoint for old build instructions.

The real Python 3.13 build configuration lives in setup.py/pyproject.toml.
"""

from pathlib import Path
from runpy import run_path


run_path(str(Path(__file__).with_name("setup.py")), run_name="__main__")
