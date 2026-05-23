import numpy as np
import pytest

import pygft


def test_real_partitions_match_matlab_core_examples():
    np.testing.assert_array_equal(pygft.real_partitions(8), np.array([1, 2, 4, 6, 8], dtype=np.int32))
    np.testing.assert_array_equal(
        pygft.real_partitions(16),
        np.array([1, 2, 4, 8, 11, 14, 16], dtype=np.int32),
    )


def test_partitions_return_int32_positive_boundaries():
    pars = pygft.partitions(16)
    assert pars.dtype == np.int32
    assert pars.ndim == 1
    assert np.all(pars > 0)


@pytest.mark.parametrize("fn", [pygft.gft1d, pygft.gft1d_real])
def test_1d_transforms_return_finite_complex_vectors(fn):
    signal = np.linspace(0.0, 1.0, 16)
    out = fn(signal)
    assert out.shape == signal.shape
    assert out.dtype == np.complex128
    assert np.all(np.isfinite(out))


def test_2d_transform_returns_finite_complex_matrix():
    image = np.arange(64, dtype=float).reshape(8, 8)
    out = pygft.gft2d(image)
    assert out.shape == image.shape
    assert out.dtype == np.complex128
    assert np.all(np.isfinite(out))


def test_interpolate_nn_returns_finite_square_output():
    signal = pygft.gft1d(np.linspace(0.0, 1.0, 16))
    out = pygft.interpolate_nn(signal, 8)
    assert out.shape == (8, 8)
    assert out.dtype == np.complex128
    assert np.all(np.isfinite(out))


def test_interpolate_nn_uses_frequency_time_layout():
    signal = np.zeros(32, dtype=np.complex128)
    signal[16] = 1.0
    transform = pygft.gft1d(signal)
    spec = np.abs(pygft.interpolate_nn(transform))
    positive_freqs = spec[:16, :]
    assert positive_freqs.shape == (16, 32)
    max_time_bin = int(np.argmax(positive_freqs.max(axis=0)))
    assert abs(max_time_bin - 16) <= 2


def test_invalid_window_type_is_rejected():
    with pytest.raises(ValueError, match="window_type"):
        pygft.gft1d(np.ones(8), window_type="hann")


def test_interpolate_rejects_upsampling_that_core_cannot_represent():
    with pytest.raises(ValueError, match="less than or equal"):
        pygft.interpolate_nn(np.ones(8, dtype=np.complex128), 16)
