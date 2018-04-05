"""
Useful function/method decorators
"""
from functools import wraps

__all__ = ['needs_dataset']


def needs_dataset(*names, default=None):
    """
    Decorator for skipping methods when the needed atomic data is not available
    """
    non_ion_datasets = ['abundance', 'ip', 'ioneq']
    names = [f'_{n}' if n not in non_ion_datasets else f'{n}' for n in names]

    def decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            if any([args[0].__getattribute__(n) is None for n in names]):
                return default
            else:
                return func(*args, **kwargs)
        return func_wrapper
    return decorator
