"""A common Cartesian projection frame for the structure and diffraction views."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _unit(vector):
    vector = np.asarray(vector, dtype=float)
    if vector.shape != (3,) or not np.isfinite(vector).all():
        raise ValueError("Enter three finite direction components.")
    scale = np.abs(vector).max()
    if scale == 0:
        raise ValueError("The beam direction must have a nonzero component.")
    scaled = vector / scale
    return scaled / np.linalg.norm(scaled)


@dataclass(frozen=True)
class Orientation:
    """Screen axes and incident beam as a right-handed Cartesian frame."""

    right: np.ndarray
    up: np.ndarray
    beam: np.ndarray

    @classmethod
    def from_cartesian(cls, vector):
        beam = _unit(vector)
        # Project the least parallel Cartesian axis onto the screen. Prefer
        # x in a tie so looking along +z places +x to the right and +y up.
        axis = np.eye(3)[np.argmin(np.abs(beam))]
        right = _unit(axis - axis.dot(beam) * beam)
        return cls(right, np.cross(beam, right), beam)

    @classmethod
    def from_zone(cls, cell, uvw):
        return cls.from_cartesian(_unit(uvw) @ np.asarray(cell))

    def orbit(self, yaw_deg, pitch_deg):
        """Yaw the beam toward screen right, then pitch it toward screen up."""
        yaw, pitch = np.deg2rad([yaw_deg, pitch_deg])
        beam = np.cos(yaw) * self.beam + np.sin(yaw) * self.right
        right = _unit(np.cross(self.up, beam))
        beam = _unit(np.cos(pitch) * beam + np.sin(pitch) * self.up)
        return Orientation(right, _unit(np.cross(beam, right)), beam)

    def roll(self, degrees):
        """Rotate the screen-right axis towards screen up about the beam."""
        angle = np.deg2rad(degrees)
        right = _unit(np.cos(angle) * self.right + np.sin(angle) * self.up)
        return Orientation(right, _unit(np.cross(self.beam, right)), self.beam)

    def zone(self, cell):
        """Return direct-lattice direction ratios with largest magnitude one."""
        uvw = np.linalg.solve(np.asarray(cell).T, self.beam)
        uvw /= np.max(np.abs(uvw))
        uvw[np.abs(uvw) < 1e-12] = 0.0
        return uvw

    def label(self, cell):
        uvw = self.zone(cell)
        for scale in range(1, 25):
            indices = np.rint(scale * uvw)
            if np.max(np.abs(indices - scale * uvw)) < 1e-7:
                return "[" + " ".join(str(int(i)) for i in indices) + "]"
        return "≈ [" + " ".join(f"{i:.3f}" for i in uvw) + "]"
