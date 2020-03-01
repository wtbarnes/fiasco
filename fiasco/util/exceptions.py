"""
Custom exceptions
"""


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
