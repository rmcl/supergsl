from typing import List, Optional
from Bio.Seq import Seq

from supergsl.core.types import SuperGSLType
from supergsl.core.types.part import Part

from pyDOE import fullfact

class AssemblyFactor:
    def __init__(self, type, levels):
        self.type = type
        self.levels = levels

    @property
    def levels(self):
        return self.levels

    def get_level(self, index):
        return self.levels[index]


class AssemblyDeclaration(SuperGSLType):
    """A declaration of a desired assembly.

    Factors -
        * Part groups
        * Assembly Options

    Levels
        Each factor has a set of possible values currenly must be discrete

    """
    def __init__(self, label : Optional[str], parts : List[Part]):
        self.label : Optional[str] = label
        self.factors = self.build_factors_from_parts(parts)

    def build_factors_from_parts(self, parts : List[Part]):
        return [
            AssemblyFactor('Part', [part])
            for part in parts
        ]

    @property
    def num_designs(self):
        """Return the number of possible designs in the assembly declaration."""
        count = 1
        for factor in self.factors:
            count *= len(factor.levels)
        return count

    def get_full_factorial_designs(self):
        """Return full-factorial iterator of the assembly designs."""

        if self.num_designs > 500:
            raise Exception('AssemblyDeclaration will generate %d designs.' % self.num_designs)

        designs = fullfact([
            len(factor.levels)
            for factor in self.factors
        ])

        for design in designs:
            yield [
                self.factors[factor_index].get_level(level_index)
                for factor_index, level_index in enumerate(design)
            ]



class AssemblyResultSet(SuperGSLType):
    """The product of assembler.


    ????
    """




class Assembly(SuperGSLType):
    """Store an assembled construct."""

    def __init__(self, identifier : str, sequence : Seq, parts : List[Part]):
        self.identifier = identifier
        self.sequence = sequence
        self.parts = parts

    def get_identifier(self):
        return self.identifier

    def get_sequence(self) -> Seq:
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
