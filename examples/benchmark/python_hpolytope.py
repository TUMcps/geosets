from __future__ import annotations

import numpy as np
from scipy.optimize import linprog


class PythonHPolytope:
    def __init__(self, *, A: np.ndarray, b: np.ndarray, validate: bool = True) -> None:
        super().__init__()
        self._tol = 1e-8

        if validate:
            if A.ndim != 2:
                raise ValueError("Matrix A must be 2D.")
            if b.ndim != 1:
                raise ValueError("Vector b must be 1D.")

            if A.shape[0] != b.shape[0]:
                raise ValueError("Matrix A must have the same number of rows as vector b.")

        self.A = A.copy()
        self.b = b.copy()
        self._center: np.ndarray | None = None

    def __str__(self) -> str:
        return f"HPolytope(A={self.A}, b={self.b})"

    @staticmethod
    def from_random(dimension: int, num_constraints: int) -> PythonHPolytope:  # pyright: ignore[reportIncompatibleMethodOverride]
        box_A = np.vstack([np.eye(dimension), -np.eye(dimension)])
        box_b = np.ones(2 * dimension)

        # Generate random hyperplanes
        random_A = np.random.standard_normal(size=(num_constraints, dimension))
        random_A /= np.linalg.norm(random_A, axis=1)[:, np.newaxis]

        interior_point = np.random.uniform(-0.8, 0.8, size=dimension)

        offsets = np.random.uniform(0.1, 1.0, size=num_constraints)
        random_b = np.dot(random_A, interior_point) + offsets

        A = np.vstack([box_A, random_A])
        b = np.hstack([box_b, random_b])

        return PythonHPolytope(A=A, b=b, validate=False)

    @staticmethod
    def from_unit_box(dimension: int) -> PythonHPolytope:
        return PythonHPolytope(
            A=np.concatenate([np.eye(dimension), -np.eye(dimension)]),
            b=np.ones(dimension * 2),
            validate=False,
        )

    def copy(self) -> PythonHPolytope:
        return PythonHPolytope(A=self.A.copy(), b=self.b.copy(), validate=False)

    def dim(self) -> int:
        return self.A.shape[1]

    def empty(self) -> bool:
        c = self.b

        # constraints
        A_eq = self.A.T
        b_eq = np.zeros(self.dim())

        # solve linear program (bounds default (0, Inf) which is required here)
        res = linprog(c, A_eq=A_eq, b_eq=b_eq)

        # cannot be infeasible since res.x = 0 fulfills any system of equalities A_eq * res.x = 0
        if res.status == 3:
            return True
        return res.fun < 0

    def center(self) -> np.ndarray:
        """Computation of the Chebyshev center of an HPolyhedron HP via linear programming.
        Defined as the center with the ball of largest radius contained in HP.

        Raises:
            EmptySetError: Set is empty.
            UnboundedSetError: Set is unbounded. LP could not converge.

        Returns:
            np.ndarray: Chebyshev center.
        """
        if self._center is not None:
            return self._center

        # objective function
        c = np.hstack((-1, np.zeros(self.dim())))

        # inequality constraints
        A_ub = np.vstack(
            (
                np.hstack((-1, np.zeros(self.dim()))),
                np.hstack(
                    (
                        np.reshape(
                            np.linalg.norm(self.A, axis=1, ord=2), (self.number_constraints(), 1)
                        ),
                        self.A,
                    )
                ),
            )
        )
        b_ub = np.hstack((0, self.b))

        # solve linear program
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(None, None))

        # check empty and unbounded cases
        if res.status == 2 or res.status == 3:  # infeasible -> empty
            raise ValueError("center failed")

        self._center = res.x[1:]
        return self._center  # type: ignore

    def number_constraints(self) -> int:
        """Returns the number of constraints in the constraint matrix.

        Returns:
            int: Number of constraints.
        """
        return self.b.size

    def translate_(self, vector: np.ndarray):
        if vector.ndim != 1:
            raise ValueError("Vector must be 1D.")
        if vector.size != self.dim():
            raise ValueError("Vector size must match the dimension of the polytope.")

        # self.b = self.b + np.linear_transform(self.A, vector)
        self.b = self.b + self.A @ vector
        self._center = None  # reset cached center

    def linear_transform_(self, matrix: np.ndarray) -> None:
        m, n = matrix.shape
        if m == n and np.linalg.matrix_rank(matrix) == n:
            # simple formula for square and invertible matrices
            # return HPolytope(A = np.linear_transform(self.A, np.linalg.inv(matrix)), b = self.b, validate = False)
            # self.A = np.linear_transform(self.A, np.linalg.inv(matrix))
            self.A = self.A @ np.linalg.inv(matrix)
            # self.b = np.linear_transform(self.b, np.linalg.inv(matrix))
        elif m > n:
            # projections to higher-dimensional space not supported
            raise NotImplementedError

    def contains_point(self, point: np.ndarray) -> bool:
        return bool(np.all(np.linear_transform(self.A, point) <= self.b))

    def boundary_point(self, direction: np.ndarray) -> np.ndarray:
        zero_centered_HP = PythonHPolytope(
            A=self.A, b=self.b - self.A @ self.center(), validate=False
        )

        if zero_centered_HP.degenerate():
            raise NotImplementedError

        # LP formulation for boundary point computation
        # min_{x,l} -l
        # s.t.      A x <= b
        #           - x + l*dir = 0

        # read out dimension
        n, h = zero_centered_HP.dim(), zero_centered_HP.number_constraints()

        # objective function
        c = np.hstack((np.zeros(n), -1.0))

        # constraints
        A_ub = np.hstack((zero_centered_HP.A, np.zeros((h, 1))))
        b_ub = zero_centered_HP.b
        A_eq = np.hstack((-np.eye(n), np.reshape(direction, (n, 1))))
        b_eq = np.zeros(n)

        # solve linear program
        res = linprog(c, A_ub, b_ub, A_eq, b_eq, bounds=(None, None))

        if res.status == 3:
            raise ValueError("Polytope is unbounded in given direction.")

        return -res.fun * direction + self.center()

    # def boundary_point_from_point_direction(
    #     self, point: np.ndarray, direction: np.ndarray
    # ) -> np.ndarray:
    #     if not self.contains_point(point):
    #         raise ValueError('The starting point is outside the polytope.')

    #     # --- Calculation ---
    #     # For each hyperplane constraint A_i * x <= b_i, we want to find the parameter 't' such that:
    #     # A_i * (point + t * direction) = b_i
    #     # t * (A_i @ direction) = b_i - (A_i @ point)
    #     # t = (b_i - A_i @ point) / (A_i @ direction)

    #     numerators = self.b - (self.A @ point)
    #     denominators = self.A @ direction

    #     valid_indices = np.where(denominators > self._tol)[0]

    #     # if valid_indices.size == 0:
    #     #     # If there are no hyperplanes that the ray moves towards, the polytope is unbounded in this direction.
    #     #     raise UnboundedSetError()

    #     # Calculate 't' for all intersecting hyperplanes
    #     t_values = numerators[valid_indices] / denominators[valid_indices]

    #     min_t = np.min(t_values)
    #     boundary_point = point + min_t * direction

    #     return boundary_point

    def line_intersection(self, point: np.ndarray, direction: np.ndarray) -> np.ndarray:
        # Vectorized computation of constraint evaluations
        # For line: x = point + t * direction
        # Constraint: A @ x <= b becomes: A @ (point + t * direction) <= b
        # Rearranged: t * (A @ direction) <= b - A @ point

        d = self.A @ direction  # (h,) - constraint coefficients for t
        c = self.b - self.A @ point  # (h,) - constraint constants

        # Handle constraints where direction is parallel to constraint hyperplane
        parallel_mask = np.abs(d) < self._tol
        if np.any(c[parallel_mask] < -self._tol):
            # Line violates parallel constraints - no intersection
            return np.empty((0, len(point)))

        # Active constraints (not parallel)
        active_mask = ~parallel_mask
        d_active = d[active_mask]
        c_active = c[active_mask]

        # Split into upper and lower bounds
        upper_mask = d_active > self._tol  # t <= c_i / d_i
        lower_mask = d_active < -self._tol  # t >= c_i / d_i

        # Compute bounds (no restriction on t >= 0 now)
        t_upper = np.inf
        t_lower = -np.inf

        if np.any(upper_mask):
            t_upper = np.min(c_active[upper_mask] / d_active[upper_mask])

        if np.any(lower_mask):
            t_lower = np.max(c_active[lower_mask] / d_active[lower_mask])

        # Check intersection validity
        if t_lower > t_upper + self._tol:
            # No intersection
            return np.empty((0, len(point)))
        elif abs(t_lower - t_upper) < self._tol:
            # Single intersection point (tangent)
            t_intersect = (t_lower + t_upper) / 2
            intersection_point = point + t_intersect * direction
            return intersection_point.reshape(1, -1)
        else:
            # Two intersection points
            intersections = np.array([point + t_upper * direction, point + t_lower * direction])
            return intersections

    def tangent_hyperplane(self, direction: np.ndarray) -> tuple[np.ndarray, float]:
        c = self.center()

        # 1. Calculate direction vector
        norm_d = np.linalg.norm(direction)
        assert norm_d > self._tol

        denominators = self.A @ direction
        numerators = self.b - (self.A @ c)

        # Use a small tolerance for floating point comparison > 0
        valid_indices = np.where(denominators > self._tol)[0]

        t_values = numerators[valid_indices] / denominators[valid_indices]
        min_t_local_idx = np.argmin(t_values)
        intersecting_halfspace_index = valid_indices[min_t_local_idx]

        return self.A[intersecting_halfspace_index], float(self.b[intersecting_halfspace_index])

    def degenerate(self, *, rtol: float = 1e-5, atol: float = 1e-8) -> bool:
        """Check if an HPolyhedron is degenerate.

        Args:
            rtol (float, optional): Relative tolerance. Defaults to 1e-5.
            atol (float, optional): Absolute tolerance. Defaults to 1e-8.

        Returns:
            bool: Degeneracy.
        """
        # compute Chebyshev center and check if it fulfills any inequality with equality
        try:
            c = self.center()
        except ValueError:
            # we consider empty sets to be degenerate
            return True

        return bool(
            np.any(np.isclose(self.b - np.linear_transform(self.A, c), 0.0, rtol=rtol, atol=atol))
        )

    def project_affine_hull(self) -> tuple:
        """Projects an HPolyhedron onto its own affine hull.
        For a degenerate HPolyhedron, the resulting HPolyhedron is of lower dimension, but non-degenerate.

        Returns:
            tuple: Projected HPolyhedron, projection matrix, center of new coordinate system in old coordinate system.
        """
        # compute basis of affine hull
        n, c = self.dim(), self.center()
        HP_shifted = self - c
        M_proj, r = HP_shifted.basis_affine_hull()
        # early exit if basis of affine hull is n-dimensional identity
        if r == n:
            return (self, np.eye(n), np.zeros(n))

        # map polyhedron onto lower dimensional space
        HP_proj = HP_shifted.linear_transform(M_proj.T)
        # remove dimensions from r to n
        A_new = HP_proj.A[:, :r]
        b_new = HP_proj.b
        # remove potential all-zero constraints
        non_flat = np.invert(np.all(np.isclose(A_new, 0.0), axis=1))
        A_new = A_new[non_flat, :]
        b_new = b_new[non_flat]
        # init resulting polyhedron
        HP_proj = PythonHPolytope(A=A_new, b=b_new)

        return (HP_proj, M_proj, c)

    def support_function(self, direction: np.ndarray) -> tuple[float, np.ndarray | None]:
        res = linprog(-direction, A_ub=self.A, b_ub=self.b, bounds=(None, None))
        if res.status == 3:  # unbounded
            return (np.inf, None)
        elif res.status == 2:  # infeasible
            return (-np.inf, None)

        return (-res.fun, res.x)

    def affine_map_to_interval_(
        self, orig_interval: np.ndarray, target_interval: np.ndarray
    ) -> None:
        """
        Affinely maps the polytope from orig_box = [[a1, a2, ...], [b1, b2, ...]]
        to target_box = [[c1, c2, ...], [d1, d2, ...]]
        """
        a, b = orig_interval
        c, d = target_interval
        a = np.asarray(a)
        b = np.asarray(b)
        c = np.asarray(c)
        d = np.asarray(d)
        assert a.shape == b.shape == c.shape == d.shape

        # Compute scaling and translation
        S = np.diag((d - c) / (b - a))
        t = c - np.diag(S) * a

        # Compute S^{-1}
        S_inv = np.diag(1.0 / np.diag(S))

        # Transform A and b
        A_new = self.A @ S_inv
        b_new = self.b + self.A @ S_inv @ t
        self.A, self.b = A_new, b_new
