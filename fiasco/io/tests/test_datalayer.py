"""
Tests for the HDF5 datalayer
"""
import pytest

from fiasco.io import DataIndexer
from fiasco.util.exceptions import MissingDatabaseError


def test_missingdatabase():
    with pytest.raises(MissingDatabaseError):
        DataIndexer.create_indexer('foo/bar.h5', '/')
