"""
Custom exceptions for fiasco.
"""

__all__ = [
    'MissingIonError',
    'MissingDatabaseError',
    'MissingASCIIFileError',
    'MissingDatasetException',
    'UnsupportedVersionError'
]

class MissingIonError(Exception):
    """
    An error to raise if an ion cannot be found in the database
    """
    pass


class MissingDatabaseError(Exception):
    """
    An error to raise when the database file cannot be found.
    """


class MissingASCIIFileError(Exception):
    """
    An error to raise when one of the CHIANTI ASCII files cannot
    be found.
    """


class MissingDatasetException(Exception):
    """
    An error to raise when a dataset file is missing.
    """


class UnsupportedVersionError(Exception):
    """
    An error to raise when an unsupported version of the CHIANTI database
    is passed to fiasco
    """
