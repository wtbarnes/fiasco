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
    names = [f'_{n}' for n in names]

    def decorator(func):
        """
        func is a method of an object that requires a dataset, e.g Ion.
        """
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            obj = args[0]
            for n in names:
                try:
                    _ = obj.__getattribute__(n)
                except KeyError:
                    raise MissingDatasetException(
                        f"{n} dataset missing for {getattr(obj, 'ion_name', obj)}."
                    )

            return func(*args, **kwargs)
        return func_wrapper
    return decorator
