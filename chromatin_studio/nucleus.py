"""Whole-nucleus envelope and genomic tethering elements.

The envelope is an optional soft ellipsoidal wall with its own centre and
semi-axes, independent of the reference confinement sphere that bounds the
modelled chromatin territory. Tethering elements are genomic intervals whose
beads can attach, one bead at a time, to an anchor point on the envelope
(lamina-like) or on the surface of an enabled nuclear-body field. An attached
bead feels a harmonic spring to its fixed anchor; attachment and detachment are
stochastic hazards evaluated in the reaction step.

These are coarse phenomenological rules: no lamina proteins, pore complexes,
membrane mechanics or anchor diffusion are represented, and the hazards are
model inputs rather than measured binding constants.
"""
import math

import numpy as np


class NucleusEnvelope:
    """Soft ellipsoid wall.  U = k/2 |u|^2 (1 - 1/rho)^2 for rho > 1.

    u = x - centre and rho = |u / a| (elementwise semi-axes a). For a sphere
    this reduces exactly to k/2 (|u| - R)^2, the same law as the reference
    confinement wall. Forces are the exact negative gradient of U.
    """

    def __init__(self, config):
        self.center = np.asarray(config.get("center", [0.0, 0.0, 0.0]), dtype=float)
        self.radii = np.asarray(config["radii"], dtype=float)
        self.stiffness = float(config.get("envelope_stiffness", 30.0))
        if self.radii.shape != (3,) or np.any(self.radii <= 0) or not np.all(np.isfinite(self.radii)):
            raise ValueError("Nucleus radii must be three positive finite semi-axes.")
        if not math.isfinite(self.stiffness) or self.stiffness < 0:
            raise ValueError("Envelope stiffness must be finite and nonnegative.")

    def rho(self, points):
        u = np.atleast_2d(points) - self.center
        return np.sqrt(np.sum((u / self.radii) ** 2, axis=1))

    def surface_point(self, points):
        """Radial projection of each point onto the envelope (along u from the centre)."""
        u = np.atleast_2d(points) - self.center
        rho = np.maximum(self.rho(points), 1e-12)
        return self.center + u / rho[:, None]

    def signed_distance(self, points):
        """Distance to the radial surface point; negative inside, positive outside."""
        u = np.atleast_2d(points) - self.center
        rho = np.maximum(self.rho(points), 1e-12)
        return np.linalg.norm(u, axis=1) * (1 - 1 / rho)

    def force_energy(self, points):
        points = np.atleast_2d(points)
        force = np.zeros_like(points)
        if self.stiffness == 0 or not len(points):
            return force, 0.0
        u = points - self.center
        rho = self.rho(points)
        out = rho > 1
        if not np.any(out):
            return force, 0.0
        uo, ro = u[out], rho[out]
        u2 = np.einsum("ij,ij->i", uo, uo)
        g = 1 - 1 / ro
        grad_rho = (uo / self.radii ** 2) / ro[:, None]
        gradient = (self.stiffness * u2 * g / ro ** 2)[:, None] * grad_rho + (self.stiffness * g ** 2)[:, None] * uo
        force[out] = -gradient
        return force, float(0.5 * self.stiffness * np.sum(u2 * g ** 2))

    def curvature_bound(self):
        """Conservative bound on the wall Hessian for integration substeps."""
        return self.stiffness * float(np.max(self.radii) / np.min(self.radii)) ** 2

    def state(self, points):
        rho = self.rho(points) if len(points) else np.zeros(0)
        _, energy = self.force_energy(points)
        return {"kind": "soft_harmonic_ellipsoid", "center": self.center.tolist(), "radii": self.radii.tolist(),
                "stiffness": self.stiffness, "outside_count": int(np.count_nonzero(rho > 1)),
                "max_scaled_radius": float(np.max(rho)) if len(rho) else 0.0, "energy": energy}


def body_surface_point(body, point):
    centre = np.asarray(body["center"], dtype=float)
    delta = np.asarray(point, dtype=float) - centre
    norm = float(np.linalg.norm(delta))
    direction = delta / norm if norm > 1e-12 else np.array([1.0, 0.0, 0.0])
    return centre + float(body["radius"]) * direction


def body_surface_distance(body, points):
    centre = np.asarray(body["center"], dtype=float)
    return np.linalg.norm(np.atleast_2d(points) - centre, axis=1) - float(body["radius"])
