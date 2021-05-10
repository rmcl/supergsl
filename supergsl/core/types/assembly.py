from typing import List, Optional
from Bio.Seq import Seq

from supergsl.core.types import SuperGSLType
from supergsl.core.types.part import Part


class AssemblyDeclaration(SuperGSLType):
    """A declaration of a desired assembly."""

    def __init__(self, label : Optional[str], parts : List[Part]):
        self.label : Optional[str] = label
        self.parts : List[Part] = parts

    def get_parts(self):
        return self.parts

class Assembly(SuperGSLType):
    """Store an assembled construct."""

    def __init__(self, identifier : str, sequence : Seq, parts : List[Part]):
        self.identifier = identifier
        self.sequence = sequence
        self.parts = parts

    def get_identifier(self):
        return self.identifier

    def get_sequence(self):
        """Return the complete sequence of the construct."""
        return self.sequence

    def get_required_parts(self):
        """Return a list of parts required to construct this assembly."""
        return self.parts

    def get_part(self):
        """Retrieve a Part corresponding to this construct."""
        raise NotImplementedError('Subclass to implement.')


class AssemblyList(SuperGSLType):
    """A collection of Assemblies."""
    def __init__(self, assemblies : List[Assembly]):
        self.assemblies = assemblies

    def __iter__(self):
        return iter(self.assemblies)
