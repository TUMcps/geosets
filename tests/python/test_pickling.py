import pickle

import numpy as np
import pytest
from definitions import GEOSET_CLASSES, get_dim_range, vertices_to_set

from geosets_py import NotImplemented


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_pickle_roundtrip_preserves_geometry(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        translated = s.translate(np.full(dim, 0.25))

        serialized = pickle.dumps(translated)
        restored = pickle.loads(serialized)

        assert isinstance(restored, SetType)
        assert restored.dim() == translated.dim()
        assert np.allclose(restored.center(), translated.center())

        try:
            assert vertices_to_set(restored.to_vertices()) == vertices_to_set(
                translated.to_vertices()
            )
        except NotImplemented:
            pytest.skip()
