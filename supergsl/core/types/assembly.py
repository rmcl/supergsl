"""Define SuperGSL Types related to the construction of new genetic assemblies"""
from typing import List, Optional, Union
from Bio.Seq import Seq
from pyDOE2 import fullfact

from supergsl.core.types import SuperGSLType
from supergsl.core.types.builtin import Collection
from supergsl.core.types.part import Part

# pylint: disable=E1136


class AssemblyFactor:
    def __init__(self, factor_type, levels):
        self._factor_type = factor_type
        self._levels = levels

    @property
    def factor_type(self):
        return self._factor_type

    @property
    def levels(self):
        return self._levels

    def __getitem__(self, index: int):
        return self._levels[index]


AssemblyDeclarationItems = Union[Part, Collection]

class AssemblyDeclaration(SuperGSLType):
    """A declaration of a desired assembly.

    Factors -
        * Part groups
        * Assembly Options

    Levels
        Each factor has a set of possible values currenly must be discrete

    """
    def __init__(
        self,
        label : Optional[str],
        items : List[AssemblyDeclarationItems]
    ):
        self._label : Optional[str] = label
        self._factors = self._build_factors_from_parts(items)

    def _build_factors_from_parts(self, items : List[AssemblyDeclarationItems]):
        factors : List[AssemblyFactor] = []

        for item in items:
            if isinstance(item, Collection):
                levels = list(item)
            else:
                levels = [item]

            factors.append(AssemblyFactor('Part', levels))
        return factors

    @property
    def label(self):
        """Return the label of this assembly declaration."""
        return self._label

    @property
    def num_designs(self):
        """Return the number of possible designs in the assembly declaration."""
        count = 1
        for factor in self._factors:
            count *= len(factor.levels)
        return count

    def get_factors(self):
        """Return all factors associated with this declaration."""
        return self._factors

    def get_levels_by_factor_type(self, factor_type : str):
        """Return the available levels for a particular factor type.

        For example, you might retrieve all parts associated with this design by
        using `declaration.get_levels_by_factor_type('Part')`.
        """
        levels = set()
        for factor in self._factors:
            if factor.factor_type != factor_type:
                continue
            levels.update(factor.levels)

        return levels

    def get_full_factorial_designs(self):
        """Return full-factorial iterator of the assembly designs."""

        if self.num_designs > 500:
            raise Exception('AssemblyDeclaration will generate %d designs.' % self.num_designs)

        designs = fullfact([
            len(factor.levels)
            for factor in self._factors
        ])

        for design in designs:
            yield [
                self._factors[factor_index][int(level_index)]
                for factor_index, level_index in enumerate(design)
            ]


class Assembly(SuperGSLType):
    """Store an assembled construct."""

    def __init__(self, identifier : str, sequence : Seq, parts : List[Part]):
        self._identifier = identifier
        self._sequence = sequence
        self._parts = parts

    @property
    def identifier(self) -> str:
        """Return the identifier of this assembly."""
        return self._identifier

    @property
    def sequence(self) -> Seq:
        """Return the complete sequence of the construct."""
        return self._sequence

    @property
    def parts(self) -> List[Part]:
        """Return a list of parts required to construct this assembly."""
        return self._parts

    def get_part(self) -> Part:
        """Retrieve a `Part` corresponding to this construct."""
        raise NotImplementedError('Subclass to implement.')


class AssemblyResultSet(SuperGSLType):
    """A collection of Assemblies."""
    def __init__(self, assemblies : List[Assembly]):
        self.assemblies = assemblies

    def __iter__(self):
        return iter(self.assemblies)

    def print(self):
        result = 'Assembly Result Set: %d assemblies\n' % len(self.assemblies)
        for assembly in self.assemblies:
            result += '    %s\n' % assembly.identifier
        return result