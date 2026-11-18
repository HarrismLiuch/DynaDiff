"""Elastic single-crystal diffraction in reciprocal cycles per ångström."""

from collections import Counter
from concurrent.futures import CancelledError
from dataclasses import dataclass
from math import prod

import gemmi
import numpy as np
from scipy import constants, fft, linalg, special


MAX_DYNAMICAL_BEAMS = 9999
_DENSE_BEAM_LIMIT = 401
# Several complex FFT arrays may coexist. Bound peak memory on a laptop.
_MAX_FFT_POINTS = 4_000_000


def _check_cancelled(cancelled):
    if cancelled is not None and cancelled():
        raise CancelledError("Simulation superseded.")


@dataclass(frozen=True)
class SimulationSettings:
    mode: str = "electron_kinematic"
    voltage_kv: float = 200.0
    wavelength_angstrom: float = 1.5406
    thickness_nm: float = 20.0
    temperature_k: float = 293.15
    g_max: float = 1.25
    max_beams: int = 101
    debye_waller: bool = True
    multiple_scattering: bool = True
    relativistic: bool = True
    ewald_curvature: bool = True


@dataclass(frozen=True)
class Pattern:
    hkl: np.ndarray
    xy: np.ndarray
    intensity: np.ndarray
    excitation_error: np.ndarray
    wavelength_angstrom: float
    notes: tuple


def electron_optics(voltage_kv, relativistic=True):
    """Return vacuum wavelength in Å and the mass multiplier for tabulated f_e."""
    if not np.isfinite(voltage_kv) or voltage_kv <= 0:
        raise ValueError("Enter a positive finite accelerating voltage.")
    energy = voltage_kv * 1000.0 * constants.e
    ratio = energy / (constants.m_e * constants.c**2) if relativistic else 0.0
    momentum = np.sqrt(2.0 * constants.m_e * energy * (1.0 + ratio / 2.0))
    return constants.h / momentum * 1e10, 1.0 + ratio


def _centering_translations(crystal, b_factors):
    """Recognize exact conventional centering in the actual weighted motif."""
    if len(crystal.frac) < 2:
        return []
    # Quantize fractional sites on a periodic grid divisible by both 2 and 3.
    grid = 6_000_000
    xyz = np.rint(crystal.frac * grid).astype(np.int64) % grid
    labels = list(zip(crystal.numbers, crystal.occupancy, b_factors))
    motif = Counter((label, *pos) for label, pos in zip(labels, xyz))
    translations = []
    # Integer sixths represent A, B, C, I and both rhombohedral settings.
    for shift in ((0, 3, 3), (3, 0, 3), (3, 3, 0), (3, 3, 3), (4, 2, 2), (2, 4, 2)):
        offset = np.asarray(shift) * (grid // 6)
        if (labels[0], *((xyz[0] + offset) % grid)) not in motif:
            continue
        moved = (xyz + offset) % grid
        if motif == Counter((label, *pos) for label, pos in zip(labels, moved)):
            translations.append(shift)
    return translations


def _beam_order(hkl, error, norm2):
    # 000 first, then near-Bragg beams, then lower frequencies. Round-off in
    # equivalent directions should not change the priority of a shell.
    return np.lexsort((hkl[:, 2], hkl[:, 1], hkl[:, 0], norm2.round(12),
                       np.abs(error).round(12), np.any(hkl != 0, axis=1)))


def _reflections(cell, basis, wavelength, settings, translations, cancelled=None):
    reciprocal = np.linalg.inv(cell).T
    # h_i = g · a_i, so these bounds also cover highly oblique cells.
    bounds = [int(np.floor(settings.g_max * n + 1e-9)) for n in np.linalg.norm(cell, axis=1)]
    shape = [2 * bound + 1 for bound in bounds]
    grid_size = prod(shape)
    if grid_size > 8_000_000:
        raise ValueError("This cell and reciprocal radius exceed the grid limit. Reduce the reciprocal radius.")
    hkl_parts, xyz_parts, error_parts, norm_parts = [], [], [], []
    candidates = 0
    k0 = 1.0 / wavelength
    xray = settings.mode == "xray_kinematic"
    dynamical = settings.mode == "electron_dynamical"
    t = 10.0 * settings.thickness_nm
    # Keep a measurable range of the slab rocking curve, including sidelobes.
    window = max(0.10, 8.0 / max(t, 1e-10))
    reciprocal_in_frame = reciprocal @ basis.T
    block = 65_536
    for start in range(0, grid_size, block):
        _check_cancelled(cancelled)
        flat = np.arange(start, min(start + block, grid_size))
        slab = np.column_stack(np.unravel_index(flat, shape)) - bounds
        keep = np.ones(len(slab), dtype=bool)
        for shift in translations:
            keep &= (slab @ shift) % 6 == 0
        hkl = slab[keep]
        xyz = hkl @ reciprocal_in_frame
        norm2 = np.einsum("ij,ij->i", xyz, xyz)
        transverse2 = xyz[:, 0]**2 + xyz[:, 1]**2
        keep = norm2 <= settings.g_max**2 * (1.0 + 1e-12)
        if xray or settings.ewald_curvature:
            keep &= transverse2 < k0**2
        if settings.ewald_curvature:
            keep &= k0 + xyz[:, 2] > 0.0
            if xray:
                kz_out = np.sqrt(np.maximum(k0**2 - transverse2, 0.0))
                # Rationalized expression avoids cancellation near the origin.
                error = -transverse2 / (k0 + kz_out) - xyz[:, 2]
            else:
                error = -xyz[:, 2] - 0.5 * wavelength * norm2
        else:
            error = -xyz[:, 2]
        if not dynamical and not xray:
            keep &= np.abs(error) <= window
        indices = np.flatnonzero(keep)
        candidates += len(indices)
        if not len(indices):
            continue
        if dynamical and len(indices) > settings.max_beams + 1:
            # The global top N can only contain local top-N candidates.
            # Keep one extra to detect split shells at the end.
            order = _beam_order(hkl[indices], error[indices], norm2[indices])
            indices = indices[order[:settings.max_beams + 1]]
        if not dynamical and candidates > 60_000:
            raise ValueError("This setting yields over 60,000 reflections. Reduce the reciprocal radius.")
        hkl_parts.append(hkl[indices])
        xyz_parts.append(xyz[indices, :2])
        error_parts.append(error[indices])
        norm_parts.append(norm2[indices])
    hkl = np.concatenate(hkl_parts)
    xy = np.concatenate(xyz_parts)
    error = np.concatenate(error_parts)
    norm2 = np.concatenate(norm_parts)
    order = _beam_order(hkl, error, norm2)
    note = None
    if dynamical and len(order) > settings.max_beams:
        count = settings.max_beams
        # Preserve whole degenerate shells whenever possible.
        excitation_rank = np.abs(error).round(12)
        radius_rank = norm2.round(12)
        last, following = order[count - 1], order[count]
        if excitation_rank[last] == excitation_rank[following] and radius_rank[last] == radius_rank[following]:
            earlier = np.flatnonzero((excitation_rank[order[:count]] == excitation_rank[last])
                                     & (radius_rank[order[:count]] == radius_rank[last]))[0]
            if earlier > 1:
                count = earlier
        note = f"Finite basis uses {count} of {candidates} candidate beams. Increase the beam limit to check convergence."
        order = order[:count]
    return hkl[order], xy[order], error[order], reciprocal, note


def _structure_factors(hkl, reciprocal, crystal, b_factors, xray, cancelled=None):
    """Sum the Fourier components for an expanded unit cell, in bounded chunks."""
    q = hkl @ reciprocal
    stol2 = np.einsum("ij,ij->i", q, q) / 4.0
    factors = np.zeros(len(hkl), dtype=np.complex128)
    symbols = np.array(crystal.symbols)
    for symbol in np.unique(symbols):
        element = gemmi.Element(str(symbol))
        coeff = element.it92 if xray else element.c4322
        if coeff is None:
            raise ValueError(f"No {'X-ray' if xray else 'electron'} form factor is tabulated for {symbol}.")
        selected = symbols == symbol
        sites = crystal.frac[selected]
        occupancy = crystal.occupancy[selected]
        b = b_factors[selected]
        # The phase and thermal arrays are (reflections, sites).
        block = max(1, 250_000 // len(sites))
        a, exponents = np.array(coeff.a), np.array(coeff.b)
        for start in range(0, len(hkl), block):
            _check_cancelled(cancelled)
            rows = slice(start, start + block)
            s2 = stol2[rows, None]
            f = np.exp(-s2 * exponents) @ a
            if xray:
                f += coeff.c
            phase = np.exp((-2j * np.pi) * (hkl[rows] @ sites.T))
            if np.any(b):
                phase *= np.exp(-s2 * b)
            factors[rows] += f * (phase @ occupancy)
    return factors


def _compact_lattice(hkl):
    """Choose integer grid coordinates with a smaller bounding box.

    Integer shears straighten tilted zone planes. Common divisors compress
    centered lattices. Every candidate retains hkl == coordinates @ transform.
    """
    coordinates = np.array(hkl, dtype=np.int64, copy=True)
    transform = np.eye(3, dtype=np.int64)
    improved = True
    while improved:
        improved = False
        for axis in range(3):
            others = [i for i in range(3) if i != axis]
            columns = coordinates[:, others]
            values = coordinates[:, axis]
            best_span = np.ptp(values)
            if best_span == 0:
                continue
            slopes = np.linalg.lstsq(columns, values, rcond=None)[0]
            m, n = np.rint(slopes).astype(np.int64)
            shears = [(0, 0)] + [(a, b) for a in range(m-1, m+2) for b in range(n-1, n+2)]
            best = None
            for a, b in shears:
                residual = values - a * columns[:, 0] - b * columns[:, 1]
                divisor = max(1, np.gcd.reduce(residual))
                span = np.ptp(residual) // divisor
                if span < best_span:
                    best_span = span
                    best = (a, b, divisor, residual // divisor)
            if best is not None:
                a, b, divisor, compact_values = best
                coordinates[:, axis] = compact_values
                transform[others[0]] += a * transform[axis]
                transform[others[1]] += b * transform[axis]
                transform[axis] *= divisor
                improved = True
    return coordinates, transform


def _fft_exit_wave(hkl, reciprocal, crystal, b_factors, scale, diagonal, thickness, cancelled):
    """Apply the same Hermitian Bloch matrix using convolution on integer HKLs."""
    coordinates, transform = _compact_lattice(hkl)
    span = np.ptp(coordinates, axis=0)
    # Each possible signed beam difference occupies a unique FFT grid index.
    shape = tuple(fft.next_fast_len(int(2 * d + 1)) for d in span)
    size = prod(shape)
    if size > _MAX_FFT_POINTS:
        raise ValueError(f"The selected beams span an FFT grid of {size:,} points "
                         f"(workspace limit {_MAX_FFT_POINTS:,}). Reduce Max. |g| or the beam limit.")
    kernel = np.zeros(shape, dtype=np.complex128)
    max_g = np.linalg.norm(hkl @ reciprocal, axis=1).max()
    shape_array = np.array(shape)
    for start in range(0, size, 32_768):
        _check_cancelled(cancelled)
        flat = np.arange(start, min(start + 32_768, size))
        delta = np.column_stack(np.unravel_index(flat, shape))
        delta = np.where(delta <= span, delta, delta - shape_array)
        # Extra cells from FFT padding, and differences outside the sphere,
        # cannot occur between retained beams. Omit them from the embedding.
        valid = np.all(np.abs(delta) <= span, axis=1)
        differences = delta[valid] @ transform
        inside = np.linalg.norm(differences @ reciprocal, axis=1) <= 2 * max_g * (1 + 1e-12)
        ids = flat[valid][inside]
        if len(ids):
            kernel.flat[ids] = scale * _structure_factors(
                differences[inside], reciprocal, crystal, b_factors, False, cancelled)
    kernel.flat[0] = 0
    # F(-h) = conj(F(h)), so the convolution has a real Fourier spectrum.
    spectrum = fft.fftn(kernel, overwrite_x=True).real.copy()
    del kernel
    grid_indices = np.ravel_multi_index(coordinates.T, shape, mode="wrap")

    # The selected-beam matrix is a principal submatrix of this circulant
    # convolution plus the excitation diagonal. Their extremal eigenvalues
    # bound the spectrum for the Chebyshev recurrence.
    low = spectrum.min() + diagonal.min()
    high = spectrum.max() + diagonal.max()
    center, radius = 0.5 * (high + low), 0.5 * (high - low)
    incident = np.zeros(len(hkl), dtype=np.complex128)
    incident[0] = 1
    if thickness * radius == 0:
        return np.exp(1j * thickness * center) * incident
    spectrum /= radius
    diagonal = (diagonal - center) / radius
    grid = np.zeros(shape, dtype=np.complex128)

    def apply(vector):
        _check_cancelled(cancelled)
        grid.fill(0)
        grid.flat[grid_indices] = vector
        frequencies = fft.fftn(grid, overwrite_x=True)
        frequencies *= spectrum
        samples = fft.ifftn(frequencies, overwrite_x=True)
        return samples.ravel()[grid_indices] + diagonal * vector

    # exp(i z x) = J0(z) + 2 sum i^n J_n(z) T_n(x), |x| <= 1.
    # The generating function bounds |J_n(z)| by exp(z*sinh(a)-n*a).
    # Sum that geometric tail to choose an absolute error below 1e-13.
    z = thickness * radius
    tolerance = 1e-13
    a = min(4.0, np.cbrt(3 * np.log(2 / tolerance) / z))
    degree = int(np.ceil((z * np.sinh(a) + np.log(2 / (-np.expm1(-a) * tolerance))) / a))
    orders = np.arange(degree + 1)
    coeff = 2 * special.jv(orders, z) * np.array([1, 1j, -1, -1j])[orders % 4]
    prev, current = incident, apply(incident)
    result = (0.5 * coeff[0]) * incident + coeff[1] * current
    for c in coeff[2:]:
        prev, current = current, 2 * apply(current) - prev
        result += c * current
    return np.exp(1j * thickness * center) * result


def simulate(crystal, frame, settings=None, *, cancelled=None):
    """Simulate a parallel beam incident along frame.beam (or row 2 of frame)."""
    settings = settings or SimulationSettings()
    _check_cancelled(cancelled)
    if settings.mode not in {"electron_kinematic", "electron_dynamical", "xray_kinematic"}:
        raise ValueError(f"Unrecognized diffraction mode: {settings.mode}.")
    if not np.isfinite(settings.g_max) or settings.g_max <= 0:
        raise ValueError("Enter a positive finite reciprocal radius in Å⁻¹.")
    if not np.isfinite(settings.thickness_nm) or settings.thickness_nm < 0:
        raise ValueError("Enter a finite, nonnegative thickness in nm.")
    if not np.isfinite(settings.temperature_k) or settings.temperature_k < 0:
        raise ValueError("Enter a finite, nonnegative temperature in K.")
    if (not isinstance(settings.max_beams, (int, np.integer))
            or not 2 <= settings.max_beams <= MAX_DYNAMICAL_BEAMS):
        raise ValueError(f"Choose a beam limit from 2 to {MAX_DYNAMICAL_BEAMS}.")
    if hasattr(frame, "beam"):
        frame = [frame.right, frame.up, frame.beam]
    basis = np.asarray(frame, dtype=float)
    if basis.shape != (3, 3) or not np.allclose(basis @ basis.T, np.eye(3), atol=1e-7):
        raise ValueError("The viewing frame must contain three orthonormal vectors.")
    xray = settings.mode == "xray_kinematic"
    if xray:
        wavelength = settings.wavelength_angstrom
        if not np.isfinite(wavelength) or wavelength <= 0:
            raise ValueError("Enter a positive finite X-ray wavelength in Å.")
        # Neutral-atom X-ray factors count electrons. r_e converts them to Å.
        scattering_scale = constants.physical_constants["classical electron radius"][0] * 1e10
    else:
        wavelength, scattering_scale = electron_optics(settings.voltage_kv, settings.relativistic)
    notes = []
    b_factors = np.zeros(len(crystal.frac))
    if settings.debye_waller:
        reference = getattr(crystal, "reference_temperature_k", 293.15)
        measured = np.isfinite(crystal.b_iso)
        # The fallback is defined at room temperature, independently of the CIF.
        b_factors = np.where(measured, np.maximum(0.0, crystal.b_iso) / reference,
                            0.5 / 293.15) * settings.temperature_k
        if not measured.all():
            notes.append(f"{(~measured).sum()} sites use estimated B = 0.5 Å² at 293 K.")
    translations = _centering_translations(crystal, b_factors)
    hkl, xy, error, reciprocal, limit_note = _reflections(
        crystal.cell, basis, wavelength, settings, translations, cancelled
    )
    if limit_note:
        notes.append(limit_note)
    scale = scattering_scale * wavelength / abs(np.linalg.det(crystal.cell))
    thickness = 10.0 * settings.thickness_nm
    if thickness == 0.0:
        intensity = np.zeros(len(hkl))
        intensity[0] = 1.0
    elif settings.mode == "electron_dynamical" and settings.multiple_scattering and len(hkl) > _DENSE_BEAM_LIMIT:
        if settings.g_max > 2 and np.linalg.norm(hkl @ reciprocal, axis=1).max() > 2:
            notes.append("Beam differences above 4 Å⁻¹ extrapolate the electron form factors.")
        amplitudes = _fft_exit_wave(hkl, reciprocal, crystal, b_factors, scale,
                                    2 * np.pi * error, thickness, cancelled)
        intensity = np.abs(amplitudes)**2
    elif settings.mode == "electron_dynamical" and settings.multiple_scattering:
        diff = hkl[:, None, :] - hkl[None, :, :]
        unique, inverse = np.unique(diff.reshape(-1, 3), axis=0, return_inverse=True)
        if settings.g_max > 2 and np.any(np.linalg.norm(unique @ reciprocal, axis=1) > 4.0 + 1e-9):
            notes.append("Some beam-pair couplings extrapolate the electron form factors beyond 4 Å⁻¹.")
        f = _structure_factors(unique, reciprocal, crystal, b_factors, False, cancelled)
        matrix = (scale * f[inverse]).reshape(len(hkl), len(hkl))
        # Discard the uniform mean-potential phase in every beam.
        np.fill_diagonal(matrix, 2.0 * np.pi * error)
        _check_cancelled(cancelled)
        values, vectors = linalg.eigh(matrix, check_finite=False)
        phases = np.exp(1j * thickness * values)
        amplitudes = vectors @ (phases * vectors[0].conj())
        intensity = np.abs(amplitudes)**2
    else:
        f = _structure_factors(hkl, reciprocal, crystal, b_factors, xray, cancelled)
        intensity = np.abs(scale * thickness * f * np.sinc(thickness * error))**2
        if xray:
            # Unpolarized Thomson polarization and exit-beam normal flux for
            # a slab at normal incidence. The grid contains only forward waves.
            sin2 = wavelength**2 * np.einsum("ij,ij->i", xy, xy)
            cos2 = np.maximum(1.0 - sin2, np.finfo(float).eps)
            intensity *= (1.0 + cos2) / (2.0 * np.sqrt(cos2))
        intensity[0] = 1.0
    return Pattern(hkl, xy, intensity, error, wavelength, tuple(notes))
