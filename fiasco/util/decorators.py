"""
Useful function/method decorators
"""
__all__ = ['has_dataset']


def has_dataset(*names, default=None):
    def decorator(func):
        def func_wrapper(*args, **kwargs):
            if any([args[0].__getattribute__('_{}'.format(n)) is None for n in names]):
                return default
            else:
                return func(*args, **kwargs)
        return func_wrapper
    return decorator
