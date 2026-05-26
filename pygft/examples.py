"""Small helpers for interactive examples.

The runnable cell-mode examples live in the top-level ``examples/`` directory.
"""

from __future__ import annotations

import numpy as np


def gaussian_window(n: int, centre: float, k: float) -> np.ndarray:
    x = np.arange(n) / n
    win = abs(k) / np.sqrt(2 * np.pi) * np.exp(-((x - 0.5) ** 2) * abs(k) ** 2 / 2)
    win = win / (np.sum(win) / n)
    return np.roll(np.roll(win, -n // 2), int(centre))
