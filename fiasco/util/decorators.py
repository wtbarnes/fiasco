"""
Useful function/method decorators.
"""
from functools import wraps

from fiasco.util.exceptions import MissingDatasetException

__all__ = ['needs_dataset']


def needs_dataset(*names):
    """
    Decorator for raising an error when the needed atomic data is not available.
    """
    non_ion_datasets = ['abundance', 'ioneq']
    names = [f'_{n}' if n not in non_ion_datasets else f'{n}' for n in names]

    def decorator(func):
        """
        func is a method of `fiasco.Ion`.
        """
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            ion = args[0]
            for n in names:
                try:
                    _ = ion.__getattribute__(n)
                except KeyError:
                    raise MissingDatasetException(f'{n} dataset missing for {ion.ion_name}.')

            return func(*args, **kwargs)
        return func_wrapper
    return decorator
