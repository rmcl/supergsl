"""Define SuperGSL Types related to the construction of new genetic assemblies"""
from typing import Generator, List, Optional, Union, Tuple, Set
from Bio.Seq import Seq
from pyDOE2 import fullfact

from supergsl.core.types import SuperGSLType
from supergsl.core.types.builtin import Collection
from supergsl.core.types.part import Part
from supergsl.core.types.position import Slice
from supergsl.core.sequence import SequenceEntry

# pylint: disable=E1136


# Assembly Level Declarations are all the types that can be used to define levels
# for a AssemblyFactor. In contrast, "Assembly Level" represents the concrete
# levels available for designs.
#
# By way of a concrete example, a "Part Collection" can be used to declare levels,
# but must be converted to its list of explicit "Parts" to be used in an
# AssemblyFactor.
AssemblyLevelDeclaration = Union[Part, Collection]
AssemblyLevel = Union[Part]

class AssemblyFactor:
    """Define a Factor and its corresponding levels for one parameter of an Assembly.

    The most concrete example of this is the part to use at an explicit position if
    the assembly. Take this example:

        let promoters = [pGAL1, pGAL3, pGAL7]
        assembly {
            uHO ; promoters ; gGENE ; dHO
        }

    The AssemblyFactor corresponding to the second position in the above assembly
    has three levels pGAL1, pGAL3, and pGAL7
    """
    def __init__(self, factor_type : str, levels : List[AssemblyLevel]):
        self._factor_type = factor_type
        self._levels : List[AssemblyLevel] = levels

    @property
    def factor_type(self) -> str:
        """Return the type of this Factor."""
        return self._factor_type

    @property
    def levels(self) -> List[AssemblyLevel]:
        """Return all possible levels for this factor."""
        return self._levels

    def __getitem__(self, index: int):
        return self._levels[index]

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
        items : List[AssemblyLevelDeclaration]
    ):
        self._label : Optional[str] = label
        self._factors : List[AssemblyFactor] = self._build_factors_from_parts(items)

    def _build_factors_from_parts(
        self,
        items : List[AssemblyLevelDeclaration]
    ) -> List[AssemblyFactor]:
        """Build a list of AssemblyFactors for the parts of this assembly declaration."""
        factors : List[AssemblyFactor] = []
        for item in items:
            levels : List[AssemblyLevel] = []
            if isinstance(item, Collection):
                levels = list(item)
            else:
                levels = [item]

            factors.append(AssemblyFactor('Part', levels))
        return factors

    @property
    def label(self) -> Optional[str]:
        """Return the label of this assembly declaration."""
        return self._label

    @property
    def num_designs(self) -> int:
        """Return the number of possible designs in the assembly declaration."""
        count = 1
        for factor in self._factors:
            count *= len(factor.levels)
        return count

    def get_factors(self) -> List[AssemblyFactor]:
        """Return all factors associated with this declaration."""
        return self._factors

    def get_levels_by_factor_type(self, factor_type : str) -> Set[AssemblyLevel]:
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

    def get_full_factorial_designs(self) -> Generator[List[AssemblyLevel], None, None]:
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

    def __init__(
        self,
        identifier : str,
        sequence_entry : SequenceEntry,
        description : Optional[str] = None
    ):
        self._identifier : str = identifier
        self._sequence_entry : SequenceEntry = sequence_entry
        self._parts : List[Tuple[Part, Slice]] = []
        self.description : Optional[str] = description

    def add_part(
        self,
        part : Part,
        target_slice : Slice
    ):
        """Add a part at a specific position in the Assembly."""

        self._parts.append(
            (part, target_slice))

    @property
    def identifier(self) -> str:
        """Return the identifier of this assembly."""
        return self._identifier

    @property
    def sequence_entry(self) -> Seq:
        """Return the sequence entry in the sequence store of the construct."""
        return self._sequence_entry

    @property
    def sequence(self) -> Seq:
        """Return the complete sequence of the construct."""
        return self._sequence_entry.sequence

    @property
    def parts_with_positions(self) -> List[Tuple[Part, Slice]]:
        """Return a list of parts and their position in this assembly's sequence."""
        return self._parts

    @property
    def parts(self) -> List[Part]:
        """Return a list of parts required to construct this assembly."""
        return [
            part_tuple[0]
            for part_tuple in self._parts
        ]

    def get_part(self) -> Part:
        """Retrieve a `Part` corresponding to this construct."""
        raise NotImplementedError('Subclass to implement.')


class AssemblyResultSet(SuperGSLType):
    """A collection of Assemblies."""
    def __init__(self, assemblies : List[Assembly]):
        self.assemblies = assemblies

    def add_assembly(self, assembly : Assembly):
        """Add an assembly to this assembly result set."""
        self.assemblies.append(assembly)

    def __iter__(self):
        """Iterating over a result set should iterate over each assembly."""
        return iter(self.assemblies)

    def print(self):
        """Return a string with details about this result set."""
        result = 'Assembly Result Set: %d assemblies\n' % len(self.assemblies)
        for assembly in self.assemblies:
            result += '    %s\n' % assembly.identifier
        return result
