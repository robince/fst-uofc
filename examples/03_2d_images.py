"""Interactive 2D GFT examples."""

# %%
import matplotlib.pyplot as plt
import numpy as np

import pygft

plt.rcParams["image.origin"] = "lower"
plt.rcParams["figure.constrained_layout.use"] = True


def radial_distance(n):
    coords = np.arange(n) - n / 2
    yy, xx = np.meshgrid(coords, coords, indexing="ij")
    return np.sqrt(xx**2 + yy**2)


def gaussian2d(n, k):
    coords = np.arange(n) / n
    win = abs(k) / np.sqrt(2 * np.pi) * np.exp(-((coords - 0.5) ** 2) * abs(k) ** 2 / 2)
    win = win / np.sum(win)
    win = np.roll(win, -n // 2)
    return np.outer(win, win)


def plot_2d_transform(image, title, window_type="gaussian"):
    n = image.shape[0]
    transform = pygft.gft2d(image, window_type=window_type)
    pars = pygft.partitions(n)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(title)
    axes[0].imshow(image.real)
    axes[0].set_title("image")
    axes[1].imshow(np.abs(transform))
    axes[1].set_title(f"|gft2d(image, {window_type!r})|")

    for p in pars[: len(pars) // 2 + 1]:
        for axis_line in (axes[1].axvline, axes[1].axhline):
            axis_line(n // 2 + p, color="tab:red", alpha=0.16)
            axis_line(n // 2 - p, color="tab:red", alpha=0.16)
    return fig


# %% 2D delta impulse
n = 128
image = np.zeros((n, n))
image[n // 2, n // 2] = 1.0

plot_2d_transform(image, "2D delta impulse: Gaussian window", "gaussian")
plot_2d_transform(image, "2D delta impulse: Box window", "box")
plt.show()

# %% Circularly symmetric cosine
r = radial_distance(n)
r = r / np.max(r)
image = np.cos(2 * np.pi * 48 * r)

plot_2d_transform(image, "Circular cosine: Gaussian window", "gaussian")
plot_2d_transform(image, "Circular cosine: Box window", "box")
plt.show()

# %% Shifted localized circular patterns
r = radial_distance(n)
r = r / np.max(r)
base = np.cos(2 * np.pi * 32 * r)
localized = base * gaussian2d(n, 8)
image = pygft.shift(localized, -n // 4, -n // 4)
image = image + pygft.shift(localized, n // 4, n // 4)

plot_2d_transform(image, "Two localized circular patterns: Gaussian window", "gaussian")
plot_2d_transform(image, "Two localized circular patterns: Box window", "box")
plt.show()
