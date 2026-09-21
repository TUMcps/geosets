from __future__ import annotations

import importlib.metadata

# Import everything
from ._gs_py_core import *

# Import of the created main nanobind module
from ._gs_py_core import (  # pyright: ignore[reportMissingModuleSource]
    GeometricSet,
    HPolytope,
    Interval,
    Polygon,
    VPolytope,
    Zonotope,
)

# Assign _plot as class method
from .plot import _plot

# Import of python code
from .python_code import print_name

GeometricSet.plot = _plot  # pyright: ignore[reportAttributeAccessIssue]

# Read the version from the pyproject.toml file
__version__ = importlib.metadata.version("geosets_py")

# Removes _gs_py_core from the path of the python objects (not necessary, but nice for the user).
GeometricSet.__module__ = "geosets_py"
Interval.__module__ = "geosets_py"
HPolytope.__module__ = "geosets_py"
VPolytope.__module__ = "geosets_py"
Polygon.__module__ = "geosets_py"

# Call to initialize randomness
random_seed(0)
