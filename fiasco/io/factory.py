"""
Factory and interface for file parser
"""
import warnings

from .sources import *
from .generic import GenericParser


class ParserFactory(type):
    """
    Metaclass for creating source classes for different CHIANTI filetypes
    """

    def __call__(cls, *args, **kwargs):
        # Use custom parser if desired
        custom_parser = None
        if 'custom_parser' in kwargs:
            custom_parser = kwargs['custom_parser']
            del kwargs['custom_parser']
        if custom_parser is not None:   
            return custom_parser(*args,**kwargs)
        # Create parser based on file extension
        filetype = args[0].split('.')[-1]
        subclass_dict = {c.filetype:c for c in all_subclasses(GenericParser) if hasattr(c,'filetype')}
        if filetype in subclass_dict:
            return subclass_dict[filetype](*args,**kwargs)
        else:
            warnings.warn('Unrecognized file extension'.format(filetype))
            return type.__call__(cls,*args,**kwargs)


def all_subclasses(cls):
    return cls.__subclasses__() + [g for s in cls.__subclasses__() for g in all_subclasses(s)]


class Parser(GenericParser, metaclass=ParserFactory):
    def __init__(self, filename, custom_parser=None, **kwargs):
        super().__init__(filename, **kwargs)
