"""Rigid fits, bounded distance calculations and geometry-only interface metrics."""

import numpy as np
from scipy.spatial.distance import cdist


def fit(moving, target):
    moving, target = np.asarray(moving), np.asarray(target)
    if moving.shape != target.shape or moving.ndim != 2 or moving.shape[1] != 3 or len(moving) < 3:
        raise ValueError("A rigid fit requires matching arrays of at least three 3D points.")
    if not np.isfinite(moving).all() or not np.isfinite(target).all():
        raise ValueError("Fit coordinates must be finite.")
    a, b = moving.mean(0), target.mean(0)
    if np.linalg.matrix_rank(moving - a) < 2 or np.linalg.matrix_rank(target - b) < 2:
        raise ValueError("Attachment atoms must not be collinear.")
    u, _, vt = np.linalg.svd((moving - a).T @ (target - b))
    correction = np.diag([1., 1., np.linalg.det(u @ vt)])
    rotation = u @ correction @ vt
    translation = b - a @ rotation
    return rotation, translation, rmsd(moving @ rotation + translation, target)


def rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((np.asarray(a) - np.asarray(b)) ** 2, axis=1))))


class Distances:
    """Float64, bounded-memory blocks; GPU never changes the scientific definition."""

    def __init__(self, device="auto"):
        self.device, self.reason, self.cp = "cpu", "CPU selected", None
        if device not in {"cpu", "auto", "cuda"}:
            raise ValueError("Device must be cpu, auto or cuda.")
        if device != "cpu":
            try:
                import cupy as cp
                if cp.cuda.runtime.getDeviceCount() < 1:
                    raise RuntimeError("No CUDA device")
                cp.sum(cp.asarray([1., 2.], dtype=cp.float64)).item()
                self.cp, self.device = cp, "cuda"
                self.reason = "CUDA distance blocks; conformer generation and fitting remain on CPU"
            except Exception as exc:
                self.reason = f"CUDA unavailable ({type(exc).__name__}); using CPU"
                if device == "cuda":
                    raise ValueError(self.reason + ". Install the GPU extra and check the NVIDIA driver.") from exc

    def blocks(self, a, b, block=128):
        # ponytail: bounded dense blocks; replace with neighbor lists for much larger assemblies.
        cp = self.cp
        b_gpu = cp.asarray(b, dtype=cp.float64) if cp is not None else None
        for start in range(0, len(a), block):
            if cp is None:
                distances = cdist(a[start:start + block], b)
            else:
                delta = cp.asarray(a[start:start + block], dtype=cp.float64)[:, None, :] - b_gpu[None, :, :]
                distances = cp.asnumpy(cp.sqrt(cp.sum(delta * delta, axis=2)))
            yield start, distances


def interface(a, b, backend, cutoff=5.0, clash=2.0):
    contacts, clashes = set(), 0
    for start, d in backend.blocks(a.xyz, b.xyz):
        clashes += int(np.count_nonzero(d < clash))
        rows, cols = np.nonzero(d <= cutoff)
        contacts.update((a.keys[start + i][:4], b.keys[j][:4]) for i, j in zip(rows, cols))
    return contacts, clashes

