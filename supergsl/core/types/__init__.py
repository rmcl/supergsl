"""Define a hierarchy of types that can be used by SuperGSL and its plugins."""

from .base import SuperGSLType
from .builtin import (
    Collection,
    NucleotideSequence,
    PrimerPair
)
