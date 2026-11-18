"""Rasterize integrated Bragg intensities on a square reciprocal-space detector."""

from functools import lru_cache

import numpy as np
from scipy.ndimage import convolve1d
from scipy.signal import fftconvolve
from scipy.special import j0, j1, ndtr


PROFILES = ("Gaussian", "Airy")
_SIGMA_PER_FWHM = 1 / np.sqrt(8 * np.log(2))
# The first positive solution of (2 J1(x) / x)**2 = 1/2.
_AIRY_HALF_MAX = 1.616339948310703


def _narrow_gaussians(x, y, weights, sigma, size):
    """Integrate narrow off-grid spots directly between detector pixel edges."""
    image = np.zeros((size, size))
    offsets = np.arange(-int(np.ceil(6 * sigma)), int(np.ceil(6 * sigma)) + 1)
    ix = np.rint(x).astype(int)[:, None] + offsets
    iy = np.rint(y).astype(int)[:, None] + offsets

    def probabilities(indices, centers):
        delta = indices - centers[:, None]
        return ndtr((delta + 0.5) / sigma) - ndtr((delta - 0.5) / sigma)

    x_weights = probabilities(ix, x) * weights[:, None]
    y_weights = probabilities(iy, y)
    for k, rows in enumerate(iy.T):
        wy = y_weights[:, k]
        valid_rows = (rows >= 0) & (rows < size)
        for m, cols in enumerate(ix.T):
            valid = valid_rows & (cols >= 0) & (cols < size)
            np.add.at(image, (rows[valid], cols[valid]), wy[valid] * x_weights[valid, m])
    return image


def _airy_kernel(width, size):
    """Integrate an isotropic Airy intensity over square pixels, including its core."""
    # The square kernel extends 16 FWHM in each direction. The full kernel
    # contains ~99% of the untruncated Airy profile's energy. A broad cropped
    # kernel represents exactly the offsets that can reach this detector.
    full_radius = int(np.ceil(16 * width))
    radius = min(full_radius, size - 1)
    axis = np.arange(radius + 1)
    edge = axis[None, :] + 0.5
    row = axis[:, None]
    # The outward flux of a radial field across a pixel edge is the angular
    # integral of the encircled energy E(r) / (2 pi). Using E = 1-J0²-J1²
    # also resolves spots smaller than one pixel without point-sampling the core.
    lower = np.arctan2(row - 0.5, edge)
    upper = np.arctan2(row + 0.5, edge)
    mid = (lower + upper) / 2
    half_span = (upper - lower) / 2
    scaled_edge = edge * (2 * _AIRY_HALF_MAX / width)
    flux = np.zeros_like(mid)
    nodes, weights = np.polynomial.legendre.leggauss(6)
    for node, weight in zip(nodes, weights):
        argument = scaled_edge / np.cos(mid + node * half_span)
        energy = 1 - j0(argument)**2 - j1(argument)**2
        flux += weight * energy
    flux *= half_span / (2 * np.pi)
    left_flux = np.concatenate((-flux[:, :1], flux[:, :-1]), axis=1)
    across = flux - left_flux
    quadrant = np.maximum(across + across.T, 0)
    half_kernel = np.concatenate((quadrant[:, :0:-1], quadrant), axis=1)
    kernel = np.concatenate((half_kernel[:0:-1], half_kernel), axis=0)
    if radius == full_radius:
        kernel /= kernel.sum()
    else:
        # Normalize by the intended finite square, even when this detector
        # only needs its central portion. On each octant its outer edge is
        # (full_radius + 1/2) / cos(theta).
        angles = (nodes + 1) * (np.pi / 8)
        argument = (full_radius + 0.5) * (2 * _AIRY_HALF_MAX / width) / np.cos(angles)
        enclosed = 1 - j0(argument)**2 - j1(argument)**2
        kernel /= np.dot(weights, enclosed) / 2
    return kernel


def detector_image(xy, intensity, g_max, *, size=513, profile="Gaussian", fwhm=0.025):
    """Return I/I₀ per detector pixel, including the direct beam.

    Pixel centers span [-g_max, +g_max] in Å⁻¹, with +qy at the top.
    FWHM is the full width through the center of one isolated spot in Å⁻¹.
    Intensities add at a shared detector position and tails can leave the image.
    """
    xy = np.asarray(xy, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    if xy.shape != (intensity.size, 2) or intensity.ndim != 1:
        raise ValueError("Expected one reciprocal-space position per intensity.")
    if not np.isfinite(xy).all() or not np.isfinite(intensity).all() or np.any(intensity < 0):
        raise ValueError("Positions and intensities must be finite and intensities nonnegative.")
    if not isinstance(size, (int, np.integer)) or size < 2:
        raise ValueError("The detector needs at least two pixels on each side.")
    if not np.isfinite(g_max) or g_max <= 0 or not np.isfinite(fwhm) or fwhm <= 0:
        raise ValueError("The detector radius and peak FWHM must be finite and positive.")
    if profile not in PROFILES:
        raise ValueError(f"Unknown detector profile: {profile}")
    keep = (intensity > 0) & (np.max(np.abs(xy), axis=1) <= g_max * (1 + 1e-12))
    if not np.any(keep):
        return np.zeros((size, size))
    scale = (size - 1) / (2 * g_max)
    x, y = np.clip((size - 1) / 2 + xy[keep] * [scale, -scale], 0, size - 1).T
    weights = intensity[keep]
    width = fwhm * scale
    sigma = width * _SIGMA_PER_FWHM
    if profile == "Gaussian" and width < 2:
        return _narrow_gaussians(x, y, weights, sigma, size)

    # Bilinear deposition retains the continuous centroid of every beam before
    # convolution, even when many beams project onto exactly the same point.
    impulses = np.zeros((size, size))
    ix, iy = x.astype(int), y.astype(int)
    fx, fy = x - ix, y - iy
    for cols, wx in ((ix, 1 - fx), (np.minimum(ix + 1, size - 1), fx)):
        for rows, wy in ((iy, 1 - fy), (np.minimum(iy + 1, size - 1), fy)):
            np.add.at(impulses, (rows, cols), weights * wx * wy)

    if profile == "Gaussian":
        radius = min(int(np.ceil(6 * sigma)), size - 1)
        centers = np.arange(-radius, radius + 1)
        kernel_1d = ndtr((centers + 0.5) / sigma) - ndtr((centers - 0.5) / sigma)
        if radius < 32:
            horizontal = convolve1d(impulses, kernel_1d, axis=1, mode="constant")
            return convolve1d(horizontal, kernel_1d, axis=0, mode="constant")
        kernel = np.outer(kernel_1d, kernel_1d)
    else:
        kernel = _airy_kernel(width, size)

    image = fftconvolve(impulses, kernel, mode="same")
    # Remove background roundoff from FFT convolution (including both signs).
    image[image < image.max() * 2e-15] = 0
    return image


def scale_image(image, *, clip_percent=0.1, logarithmic=True):
    """Map intensities to [0, 1] and return the limits in I/I₀ per pixel."""
    image = np.asarray(image, dtype=float)
    if not 0 <= clip_percent < 100:
        raise ValueError("The upper-tail clipping percentage must be in [0, 100).")
    peak = float(image.max(initial=0))
    if peak == 0:
        return np.zeros_like(image), 0.0, 0.0
    upper = float(np.percentile(image, 100 - clip_percent)) if clip_percent else peak
    if upper <= 0:
        # A very sparse pattern can have a percentile inside the empty
        # background. Scale using its positive pixels in that case.
        upper = float(np.percentile(image[image > 0], 100 - clip_percent))
    if logarithmic:
        lower = max(upper * 1e-4, np.nextafter(0.0, 1.0))
        scaled = (np.log10(np.clip(image, lower, upper)) - np.log10(lower)) / 4
    else:
        lower = 0.0
        scaled = np.clip(image, 0, upper) / upper
    return np.clip(scaled, 0, 1), lower, upper


@lru_cache(maxsize=1)
def colormap_groups():
    """Available forward maps in each provider's published gallery order."""
    from matplotlib import colormaps
    import cmcrameri.cm  # noqa: F401 -- registers scientific colour maps as cmc.*
    import cmasher  # noqa: F401 -- registers CMasher maps as cmr.*

    from .mpl_groups import matplotlib_groups
    from .scientific_groups import scientific_groups

    mpl = tuple((f"Matplotlib · {label}", names) for label, names in matplotlib_groups(colormaps))
    return mpl + scientific_groups(colormaps)


@lru_cache(maxsize=1)
def colormap_names():
    """Flatten the picker catalog for callers that only need registry names."""
    return tuple(name for _, names in colormap_groups() for name in names)


@lru_cache(maxsize=32)
def color_lut(name, reverse=False):
    """Return a 256-entry RGBA lookup table ordered from low to high signal."""
    from matplotlib import colormaps

    if name not in colormaps:
        colormap_names()
    lut = colormaps[name](np.linspace(0, 1, 256), bytes=True)
    if reverse:
        lut = lut[::-1].copy()
    lut.setflags(write=False)
    return lut
