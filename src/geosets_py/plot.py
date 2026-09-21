import matplotlib
import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np


def _plot(
    self,
    ax: matplotlib.axes.Axes | None = None,
    axis: tuple = (0, 1),
    show: bool = False,
    fill_opacity: float = 0.0,
    **kwargs,
):
    """
    Plots the convex set.
    `kwargs` are given to plt.plot.
    """
    # self._check_axis(axis)

    # project onto axes
    # projected_set = self.project(axis=axis)

    if ax is None:
        ax = plt.gca()

    # compute vertices
    vertices: np.ndarray = self.to_vertices()
    vertices = vertices[:, axis]

    # Order the vertices clockwise
    center = np.mean(vertices, axis=0)
    angles = np.arctan2(vertices[:, 1] - center[1], vertices[:, 0] - center[0])
    sort_order = np.argsort(angles)
    vertices = vertices[sort_order, :]

    vertices = np.vstack((vertices, vertices[0, :]))

    # plot
    if vertices.shape[0] == 1:
        # single point: add marker
        ax.plot(vertices[:, 0], vertices[:, 1], "o", **kwargs)
    else:
        ax.plot(vertices[:, 0], vertices[:, 1], **kwargs)
        if fill_opacity > 0.0:
            ax.fill(vertices[:, 0], vertices[:, 1], alpha=fill_opacity, **kwargs)

    if show:
        plt.show()
