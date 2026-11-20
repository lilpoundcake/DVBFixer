"""Private reduced coordinates for exact observations and FF-equilibrium X–H.

Heavy mobile atoms are Cartesian. Each mobile auxiliary H uses a stereographic
two-coordinate sphere chart about its actual heavy parent. Fixed observations
are not optimization variables. The excluded chart antipode is a documented
local-optimizer limitation, not a reason to accept a small chart gradient;
stationarity is checked independently in Cartesian tangent space.
"""
from __future__ import annotations

import math

import numpy as np
from scripts.anchor_refinement_guards import _axis, xh_equilibria

REVISION = "private-reduced-xh-stereographic-v1"
EQUALITY_TOLERANCE_NM = 1.e-12


class ReducedCoordinates:
    def __init__(self, topology, system, initial_nm):
        atoms, points, masses = _axis(topology, system, initial_nm)
        if not np.isfinite(points).all() or system.getNumConstraints():
            raise ValueError("reduced coordinates require finite unconstrained initial FF state")
        if any(system.isVirtualSite(i) for i in range(len(atoms))):
            raise ValueError("virtual sites require separate reduced-coordinate support")
        self.reference = points.copy()
        self.mobile = masses > 0
        self.heavy = np.array([a.index for a in atoms if a.element.symbol != "H" and self.mobile[a.index]], dtype=int)
        self.bonds = xh_equilibria(topology, system)
        self.h = np.array([r["hydrogen_index"] for r in self.bonds], dtype=int)
        self.parent = np.array([r["parent_index"] for r in self.bonds], dtype=int)
        self.length = np.array([r["equilibrium_nm"] for r in self.bonds])
        if not self.mobile[self.h].all():
            raise ValueError("fixed/observed hydrogens are unsupported; never move them to fit an equilibrium")
        vector = points[self.h]-points[self.parent]
        distance = np.linalg.norm(vector, axis=1)
        if (distance <= 1.e-12).any():
            raise ValueError("undefined initial X–H direction")
        self.normal = vector/distance[:, None]
        if len(self.h):
            axes = np.eye(3)[np.argmin(np.abs(self.normal), axis=1)]
            self.first = np.cross(self.normal, axes)
            self.first /= np.linalg.norm(self.first, axis=1)[:, None]
            self.second = np.cross(self.normal, self.first)
        else:
            self.first = self.second = np.empty((0, 3))
        self.heavy_size = 3*len(self.heavy)
        self.size = self.heavy_size+2*len(self.h)
        self.dof = 3*int(self.mobile.sum())-len(self.h)
        if self.size != self.dof or not self.size:
            raise ValueError("reduced coordinate DOF inventory differs")
        self.initial = np.concatenate((points[self.heavy].ravel(), np.zeros(2*len(self.h))))
        self.parent_groups = {int(p): np.flatnonzero(self.parent == p) for p in set(self.parent)}

    def _chart(self, variables):
        q = np.asarray(variables, dtype=float)
        if q.shape != (self.size,) or not np.isfinite(q).all():
            raise ValueError("invalid reduced-coordinate vector")
        uv = q[self.heavy_size:].reshape(-1, 2)/(2*self.length[:, None])
        radius2 = np.sum(uv*uv, axis=1)
        denominator = 1+radius2
        direction = (2*uv[:, :1]*self.first+2*uv[:, 1:]*self.second+(1-radius2[:, None])*self.normal)/denominator[:, None]
        if not np.isfinite(direction).all():
            raise ValueError("nonfinite stereographic chart; no silent rechart/restart")
        return q, uv, denominator, direction

    def positions(self, variables):
        q, _, _, direction = self._chart(variables)
        result = self.reference.copy()
        result[self.heavy] = q[:self.heavy_size].reshape(-1, 3)
        result[self.h] = result[self.parent]+self.length[:, None]*direction
        return result

    def pullback(self, variables, cartesian_gradient):
        """Analytic dE/dq, including translation of H with a mobile parent."""
        _, uv, denominator, direction = self._chart(variables)
        gradient = np.asarray(cartesian_gradient, dtype=float)
        if gradient.shape != self.reference.shape or not np.isfinite(gradient).all():
            raise ValueError("invalid Cartesian energy gradient")
        combined = gradient.copy()
        np.add.at(combined, self.parent, gradient[self.h])
        # q_H is in nm: uv=q_H/(2*r0), so dr_H/dq_H is half dn/duv.
        du = (self.first-uv[:, :1]*(self.normal+direction))/denominator[:, None]
        dv = (self.second-uv[:, 1:]*(self.normal+direction))/denominator[:, None]
        hydrogen = np.column_stack((np.sum(gradient[self.h]*du, axis=1), np.sum(gradient[self.h]*dv, axis=1)))
        return np.concatenate((combined[self.heavy].ravel(), hydrogen.ravel()))

    def project_forces(self, positions, forces):
        """Euclidean Cartesian nullspace projection; independent of sphere charts.

        For each heavy parent, solve (J J^T) lambda=J F, then F-J^T lambda.
        Fixed-parent blocks are identity; mobile-parent blocks are I+N N^T.
        The fixed coordinates are removed before projection, not held by springs.
        """
        points, raw = np.asarray(positions), np.asarray(forces)
        if points.shape != self.reference.shape or raw.shape != points.shape or not np.isfinite(raw).all():
            raise ValueError("invalid tangent-projection axes/forces")
        vector = points[self.h]-points[self.parent]
        lengths = np.linalg.norm(vector, axis=1)
        if not np.isfinite(lengths).all() or (lengths <= 1.e-12).any():
            raise ValueError("undefined X–H constraint normal")
        normals = vector/lengths[:, None]
        projected = raw.copy()
        projected[~self.mobile] = 0
        for parent, group in self.parent_groups.items():
            n = normals[group]
            h = self.h[group]
            parent_force = raw[parent] if self.mobile[parent] else np.zeros(3)
            rhs = np.sum(n*(raw[h]-parent_force), axis=1)
            gram = np.eye(len(group))+(n@n.T if self.mobile[parent] else 0)
            multipliers = np.linalg.solve(gram, rhs)
            corrections = multipliers[:, None]*n
            projected[h] -= corrections
            if self.mobile[parent]:
                projected[parent] += corrections.sum(axis=0)
        return projected

    def audit(self, positions, forces):
        points = np.asarray(positions)
        if points.shape != self.reference.shape or not np.isfinite(points).all():
            raise ValueError("invalid reduced-coordinate returned positions")
        residual = np.abs(np.linalg.norm(points[self.h]-points[self.parent], axis=1)-self.length)
        maximum = float(residual.max()) if residual.size else 0.
        projected = self.project_forces(points, forces)
        scale = float(np.max(np.abs(projected)))
        rms = 0. if scale == 0 else scale*math.sqrt(float(np.sum((projected/scale)**2))/self.dof)
        return {"fixed_exact": np.array_equal(points[~self.mobile], self.reference[~self.mobile]),
                "xh_max_equality_error_nm": maximum, "xh_equalities_satisfied": maximum <= EQUALITY_TOLERANCE_NM,
                "equality_numeric_tolerance_nm": EQUALITY_TOLERANCE_NM, "free_dof": self.dof,
                "projected_cartesian_dof_rms_kj_mol_nm": rms,
                "stationarity_passed": rms <= 10.,
                "force_scale_kj_mol_nm": 10., "constraint_count": len(self.h),
                "scope": "Cartesian tangent projection, not raw radial reaction forces or chart-gradient convergence"}
