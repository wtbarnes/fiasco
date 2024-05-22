"""
Factory and interface for file parser.
"""
import pathlib

from fiasco.io.generic import GenericParser
from fiasco.io.sources import *

__all__ = ['Parser']

class ParserFactory(type):
    """
    Metaclass for creating source classes for different CHIANTI filetypes
    """

    def __call__(cls, *args, **kwargs):
        # Import here to avoid circular imports
        from fiasco import log

        path = pathlib.Path(args[0])

        # Allow for standalone files
        if path.exists():
            kwargs['standalone'] = True
        # Use custom parser if desired
        custom_parser = None
        if 'custom_parser' in kwargs:
            custom_parser = kwargs['custom_parser']
            del kwargs['custom_parser']
        if custom_parser is not None:
            return custom_parser(*args, **kwargs)

        # Create parser based on file extension or name
        filetype_name = path.stem
        filetype_ext = path.suffix[1:]

        subclass_dict = {c.filetype: c for c in all_subclasses(GenericParser) if hasattr(c, 'filetype')}
        if filetype_ext in subclass_dict:
            return subclass_dict[filetype_ext](*args, **kwargs)
        elif filetype_name in subclass_dict:
            return subclass_dict[filetype_name](*args, **kwargs)
        else:
            log.debug(f'Unrecognized filename and extension {args[0]}', stacklevel=2)
            return type.__call__(cls, *args, **kwargs)


def all_subclasses(cls):
    """
    Return all subclasses of a given class
    """
    return cls.__subclasses__() + [g for s in cls.__subclasses__() for g in all_subclasses(s)]


class Parser(GenericParser, metaclass=ParserFactory):
    """
    General parser interface for all CHIANTI datatypes.

    The Parser ingests the name of a raw ASCII data file and builds an
    `astropy.table.QTable` from it. A predefined parser is created based
    on the file extension, but a custom parser can also be used.
    """
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
