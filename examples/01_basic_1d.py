"""Interactive 1D GFT examples.

Run this file cell-by-cell in VS Code, Jupyter, or any editor that supports
Python `# %%` cells.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

import pygft

plt.rcParams["image.origin"] = "lower"
plt.rcParams["figure.constrained_layout.use"] = True


def gaussian_window(n, centre, k):
    x = np.arange(n) / n
    win = abs(k) / np.sqrt(2 * np.pi) * np.exp(-((x - 0.5) ** 2) * abs(k) ** 2 / 2)
    win = win / (np.sum(win) / n)
    return np.roll(np.roll(win, -n // 2), int(centre))


def plot_transform(signal, title, window_type="gaussian"):
    n = signal.size
    x = np.arange(n) / n
    transform = pygft.gft1d(signal, window_type=window_type)
    pars = pygft.partitions(n)

    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=False)
    fig.suptitle(title)
    axes[0].plot(x, signal.real, label="real")
    axes[0].plot(x, signal.imag, "--", label="imag")
    axes[0].set_title("signal")
    axes[0].legend(loc="upper right")

    axes[1].plot(np.abs(transform))
    axes[1].set_title(f"|gft1d(signal, {window_type!r})|")
    for p in pars[: len(pars) // 2 + 1]:
        axes[1].axvline(n // 2 + p, color="tab:red", alpha=0.18)
        axes[1].axvline(n // 2 - p, color="tab:red", alpha=0.18)

    spec = pygft.interpolate_nn(transform)
    positive_freqs = np.abs(spec[: n // 2, :])
    axes[2].imshow(positive_freqs, aspect="auto")
    axes[2].set_title("nearest-neighbor interpolated spectrogram")
    axes[2].set_xlabel("time")
    axes[2].set_ylabel("frequency bin")
    return fig


# %% Delta impulse, Gaussian and box windows
n = 256
signal = np.zeros(n, dtype=np.complex128)
signal[n // 2] = 1.0

plot_transform(signal, "Delta impulse: Gaussian window", "gaussian")
plot_transform(signal, "Delta impulse: Box window", "box")
plt.show()

# %% Two impulses
signal = np.zeros(n, dtype=np.complex128)
signal[n // 4] = 1.0
signal[n // 4 + n // 2] = 1.0

plot_transform(signal, "Two impulses: Gaussian window", "gaussian")
plot_transform(signal, "Two impulses: Box window", "box")
plt.show()

# %% Gaussian-windowed sinusoid packets
x = np.arange(n) / n
signal = (
    np.sin(2 * np.pi * 24 * x) * gaussian_window(n, n // 4, 8)
    + np.sin(2 * np.pi * 96 * x) * gaussian_window(n, 3 * n // 4, 8)
    + np.sin(2 * np.pi * 80 * x) * gaussian_window(n, 3 * n // 4, 8)
)

plot_transform(signal, "Three localized sine packets: Gaussian window", "gaussian")
plot_transform(signal, "Three localized sine packets: Box window", "box")
plt.show()

# %% Chirp-like increasing-frequency sine
signal = np.sin(2 * np.pi * x * x * n / 4)

plot_transform(signal, "Increasing-frequency sine: Gaussian window", "gaussian")
plot_transform(signal, "Increasing-frequency sine: Box window", "box")
plt.show()

# %% Linear chirp with a known start and end frequency
duration = 1.0
start_hz = 8.0
end_hz = 96.0
t = np.arange(n) / n * duration
chirp_rate = (end_hz - start_hz) / duration
phase = 2 * np.pi * (start_hz * t + 0.5 * chirp_rate * t**2)
signal = np.sin(phase)

fig = plot_transform(signal, "Linear chirp: 8 Hz to 96 Hz", "gaussian")
fig.axes[0].set_xlabel("time (s)")
fig.axes[2].plot(
    np.arange(n),
    start_hz + chirp_rate * t,
    color="white",
    linewidth=1.2,
    alpha=0.85,
)
fig.axes[2].set_ylim(0, n // 2)
plt.show()
