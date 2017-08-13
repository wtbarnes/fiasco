"""
Factory and interface for file parser
"""
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
        if filetype == 'elvlc':
            return ElvlcParser(*args,**kwargs)
        elif filetype == 'wgfa':
            return WgfaParser(*args,**kwargs)
        elif filetype == 'fblvl':
            return FblvlParser(*args,**kwargs)
        elif filetype == 'scups':
            return ScupsParser(*args,**kwargs)
        elif filetype == 'easplom':
            return EasplomParser(*args,**kwargs)
        elif filetype == 'easplups':
            return EasplupsParser(*args,**kwargs)
        elif filetype == 'psplups':
            return PsplupsParser(*args,**kwargs)
        elif filetype == 'cilvl':
            return CilvlParser(*args,**kwargs)
        elif filetype == 'reclvl':
            return ReclvlParser(*args,**kwargs)
        elif filetype == 'rrparams':
            return RrparamsParser(*args,**kwargs)
        elif filetype == 'trparams':
            return TrparamsParser(*args,**kwargs)
        elif filetype == 'drparams':
            return DrparamsParser(*args,**kwargs)
        elif filetype == 'diparams':
            return DiparamsParser(*args,**kwargs)
        elif filetype == 'abund':
            return AbundParser(*args,**kwargs)
        elif filetype == 'ioneq':
            return IoneqParser(*args,**kwargs)
        else:
            return type.__call__(cls,*args,**kwargs)


class Parser(GenericParser, metaclass=ParserFactory):
    def __init__(self, filename, custom_parser=None, **kwargs):
        super().__init__(filename, **kwargs)
