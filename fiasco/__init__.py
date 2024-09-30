"""
fiasco: A Python interface to the CHIANTI atomic database
"""
from fiasco.collections import IonCollection
from fiasco.elements import Element
from fiasco.fiasco import (
    get_isoelectronic_sequence,
    list_elements,
    list_ions,
    proton_electron_ratio,
)
from fiasco.gaunt import GauntFactor
from fiasco.ions import Ion
from fiasco.levels import Level, Transitions
from fiasco.util.util import setup_paths

try:
    from fiasco.version import __version__
except ImportError:
    __version__ = "unknown"

defaults = setup_paths()


from fiasco.util.logger import _init_log

log = _init_log()


def _get_bibtex():
    import textwrap

    from itertools import compress
    from pathlib import Path

    # Set the bibtex entry to the article referenced in CITATION.rst
    citation_file = Path(__file__).parent / "CITATION.rst"

    # Explicitly specify UTF-8 encoding in case the system's default encoding is problematic
    with Path.open(citation_file, "r", encoding="utf-8") as citation:
        # Extract the first bibtex block:
        ref = citation.read().partition(".. code:: bibtex\n\n")[2]
        lines = ref.split("\n")
        # Only read the lines which are indented
        lines = list(compress(lines, [line.startswith("   ") for line in lines]))
        return textwrap.dedent("\n".join(lines))


__citation__ = __bibtex__ = _get_bibtex()

__all__ = [
    "IonCollection",
    "Element",
    "list_elements",
    "list_ions",
    "get_isoelectronic_sequence",
    "proton_electron_ratio",
    "Ion",
    "Level",
    "Transitions",
    "GauntFactor",
    "defaults",
    "log",
    "__version__",
    "__citation__",
]
