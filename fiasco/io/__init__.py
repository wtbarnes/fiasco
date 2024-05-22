"""
Parsers for raw CHIANTI atomic data.
fiasco uses these parsers to transform the raw ASCII files into an HDF5 version of the CHIANTI database.
"""
from fiasco.io.datalayer import DataIndexer, DataIndexerHDF5
from fiasco.io.factory import Parser
from fiasco.io.generic import GenericIonParser, GenericParser

__all__ = ["DataIndexer", "DataIndexerHDF5", "Parser", "GenericIonParser", "GenericParser"]
