"""
Useful function/method decorators
"""
from functools import wraps

from fiasco.util import MissingDatasetException

__all__ = ['needs_dataset']


def needs_dataset(*names):
    """
    Decorator for raising an error when the needed atomic data is not available.
    """
    non_ion_datasets = ['abundance', 'ioneq']
    names = [f'_{n}' if n not in non_ion_datasets else f'{n}' for n in names]

    def decorator(func):
        """
        func is a method of `fiasco.ion.Ion`.
        """
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            ion = args[0]
            missing = [n for n in names if ion.__getattribute__(n) is None]
            if len(missing):
                raise MissingDatasetException(f'Missing {missing} files for {ion.ion_name}.')

            return func(*args, **kwargs)
        return func_wrapper
    return decorator
