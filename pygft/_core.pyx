# cython: boundscheck=False, wraparound=False, initializedcheck=False

from libc.stdlib cimport free
from libc.string cimport memcpy

import numpy as np
cimport numpy as cnp


cnp.import_array()

ctypedef cnp.complex128_t complex128_t
ctypedef void (*windowFunction)(double *, int, int)


cdef extern from "gft.h":
    int gft_1dSizeOfPartitions(unsigned int N)
    int *gft_1dPartitions(unsigned int N)
    int *gft_1dRealPartitions(unsigned int N)
    void gaussian(double *win, int N, int freq)
    void box(double *win, int N, int freq)
    double *windows(int N, windowFunction window)
    double *windowsFromPars(int N, windowFunction window, int *pars)
    void gft_1dComplex64(double *signal, unsigned int N, double *win, int *pars, int stride)
    void gft_2dComplex64(double *image, unsigned int N, unsigned int M, windowFunction window)
    double *gft_1d_interpolateNN(double *signal, unsigned int N, unsigned int M)


cdef unsigned int _as_uint_size(object value) except? 0:
    cdef unsigned int n
    try:
        n = <unsigned int>value
    except (TypeError, ValueError):
        n = <unsigned int>np.max(np.shape(value))
    if n == 0:
        raise ValueError("size must be positive")
    return n


cdef cnp.ndarray _copy_partitions(int *ptr):
    cdef Py_ssize_t count = 0
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] result
    if ptr == NULL:
        raise MemoryError("failed to allocate partition vector")
    while ptr[count] > 0:
        count += 1
    result = np.empty(count, dtype=np.int32)
    if count:
        memcpy(<cnp.int32_t *>result.data, ptr, count * sizeof(int))
    return result


cdef cnp.ndarray[complex128_t, ndim=1, mode="c"] _as_complex_vector(object sig):
    cdef cnp.ndarray arr = np.ascontiguousarray(sig, dtype=np.complex128)
    if arr.ndim != 1:
        raise ValueError("input must be a 1D vector")
    if arr.shape[0] == 0:
        raise ValueError("input vector must not be empty")
    return arr


cdef cnp.ndarray[complex128_t, ndim=2, mode="c"] _as_complex_matrix(object image):
    cdef cnp.ndarray arr = np.ascontiguousarray(image, dtype=np.complex128)
    if arr.ndim != 2:
        raise ValueError("input must be a 2D matrix")
    if arr.shape[0] == 0 or arr.shape[1] == 0:
        raise ValueError("input matrix must not be empty")
    return arr


cdef windowFunction _window_from_name(str window_type) except NULL:
    if window_type == "gaussian":
        return gaussian
    if window_type == "box":
        return box
    raise ValueError("window_type must be 'gaussian' or 'box'")


def phase(num):
    arr = np.asarray(num, dtype=np.complex128)
    return np.arctan2(arr.imag, arr.real)


def polar(r, p):
    return complex(r * np.cos(p), r * np.sin(p))


def shift(a, nx, ny=None):
    arr = np.asarray(a)
    if arr.ndim == 1:
        return np.roll(arr, int(round(nx)))
    if arr.ndim == 2:
        if ny is None:
            raise ValueError("ny is required for 2D shifts")
        return np.roll(np.roll(arr, int(round(ny)), axis=0), int(round(nx)), axis=1)
    raise ValueError("input must be 1D or 2D")


def partitions(N):
    cdef unsigned int n = _as_uint_size(N)
    cdef int *ptr = gft_1dPartitions(n)
    try:
        return _copy_partitions(ptr)
    finally:
        free(ptr)


def real_partitions(N):
    cdef unsigned int n = _as_uint_size(N)
    cdef int *ptr = gft_1dRealPartitions(n)
    try:
        return _copy_partitions(ptr)
    finally:
        free(ptr)


def gaussian_windows(N):
    cdef unsigned int n = _as_uint_size(N)
    cdef double *win = windows(n, gaussian)
    cdef cnp.ndarray[complex128_t, ndim=1, mode="c"] result = np.empty(n, dtype=np.complex128)
    if win == NULL:
        raise MemoryError("failed to allocate windows")
    try:
        memcpy(<complex128_t *>result.data, win, n * 2 * sizeof(double))
    finally:
        free(win)
    return result.real


def gft1d(sig, window_type="gaussian"):
    cdef cnp.ndarray[complex128_t, ndim=1, mode="c"] result = _as_complex_vector(sig).copy()
    cdef unsigned int n = result.shape[0]
    cdef windowFunction window = _window_from_name(window_type)
    cdef int *pars = gft_1dPartitions(n)
    cdef double *win = windows(n, window)
    if pars == NULL or win == NULL:
        free(pars)
        free(win)
        raise MemoryError("failed to allocate transform workspace")
    try:
        gft_1dComplex64(<double *>result.data, n, win, pars, 1)
    finally:
        free(pars)
        free(win)
    return np.roll(result, n // 2)


def gft1d_real(sig):
    cdef cnp.ndarray[complex128_t, ndim=1, mode="c"] result = _as_complex_vector(sig).copy()
    cdef unsigned int n = result.shape[0]
    cdef int *pars = gft_1dRealPartitions(n)
    cdef double *win = windowsFromPars(n, gaussian, pars)
    if pars == NULL or win == NULL:
        free(pars)
        free(win)
        raise MemoryError("failed to allocate transform workspace")
    try:
        gft_1dComplex64(<double *>result.data, n, win, pars, 1)
    finally:
        free(pars)
        free(win)
    return np.roll(result, n // 2)


def gft2d(image, window_type="gaussian"):
    cdef cnp.ndarray[complex128_t, ndim=2, mode="c"] result = _as_complex_matrix(image).copy()
    cdef unsigned int n = result.shape[0]
    cdef unsigned int m = result.shape[1]
    cdef windowFunction window = _window_from_name(window_type)
    gft_2dComplex64(<double *>result.data, n, m, window)
    return np.roll(np.roll(result, n // 2, axis=0), m // 2, axis=1)


def interpolate_nn(SIG, M=None):
    cdef cnp.ndarray[complex128_t, ndim=1, mode="c"] sig = _as_complex_vector(SIG)
    cdef cnp.ndarray[complex128_t, ndim=1, mode="c"] shifted
    cdef cnp.ndarray[complex128_t, ndim=2, mode="c"] result
    cdef unsigned int n = sig.shape[0]
    cdef unsigned int m = n if M is None else _as_uint_size(M)
    cdef double *image

    if m > n:
        raise ValueError("M must be less than or equal to the input length")

    shifted = np.ascontiguousarray(np.roll(sig, -(<int>n // 2)), dtype=np.complex128)
    result = np.empty((m, m), dtype=np.complex128)
    image = gft_1d_interpolateNN(<double *>shifted.data, n, m)
    if image == NULL:
        raise MemoryError("failed to allocate interpolation output")
    try:
        memcpy(<double *>result.data, image, m * m * 2 * sizeof(double))
    finally:
        free(image)
    return result
