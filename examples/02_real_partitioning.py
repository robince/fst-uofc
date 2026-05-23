"""Interactive examples for the real-valued asymmetric partitioning path."""

# %%
import matplotlib.pyplot as plt
import numpy as np

import pygft

plt.rcParams["image.origin"] = "lower"
plt.rcParams["figure.constrained_layout.use"] = True


def plot_real_vs_complex(signal, title):
    n = signal.size
    complex_gft = pygft.gft1d(signal, window_type="gaussian")
    real_gft = pygft.gft1d_real(signal)

    fig, axes = plt.subplots(3, 1, figsize=(10, 8))
    fig.suptitle(title)

    axes[0].plot(signal)
    axes[0].set_title("real input signal")

    axes[1].plot(np.abs(complex_gft), label="gft1d")
    axes[1].plot(np.abs(real_gft), "--", label="gft1d_real")
    axes[1].set_title("transform magnitudes")
    axes[1].legend(loc="upper right")

    real_pars = pygft.real_partitions(n)
    axes[2].plot(np.abs(real_gft))
    axes[2].set_title("real-valued asymmetric partitions")
    for p in real_pars:
        axes[2].axvline((n // 2 + p) % n, color="tab:red", alpha=0.25)
    return fig


# %% Confirm the asymmetric real-valued partitions
for n in (8, 16, 32, 64):
    print(f"N={n}: {pygft.real_partitions(n).tolist()}")

# %% Delta impulse using the real-valued transform path
n = 256
signal = np.zeros(n)
signal[n // 2] = 1.0

plot_real_vs_complex(signal, "Delta impulse: complex vs real-valued GFT")
plt.show()

# %% Pure sinusoid similar to test_st.m
freq_hz = 6
n = 512
x = np.linspace(0, 2, n)
signal = np.sin(2 * np.pi * freq_hz * x)

plot_real_vs_complex(signal, "Pure sinusoid: complex vs real-valued GFT")
plt.show()

# %% Real-partitioned spectrogram
transform = pygft.gft1d_real(signal)
positive_freqs = pygft.interpolate_real(transform)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].imshow(np.abs(positive_freqs), aspect="auto")
axes[0].set_title("real-partition magnitude")
axes[0].set_xlabel("time")
axes[0].set_ylabel("frequency bin")
axes[1].imshow(np.angle(positive_freqs), aspect="auto", cmap="twilight")
axes[1].set_title("real-partition phase")
axes[1].set_xlabel("time")
plt.show()
