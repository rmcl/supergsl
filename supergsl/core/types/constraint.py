"""Define the Constraint base class for imposing design constraints on Assemblies."""
from typing import List
from supergsl.core.types.assembly import Assembly, AssemblyLevel


class Constraint:
    """Define a constraint on an AssemblyDeclaration."""

    is_definition_constraint = False
    is_assembly_constraint = False

    def evaluate_definition(self, design_definition : List[AssemblyLevel]) -> bool:
        """Evaluate the constraint based on the ordered list of parts."""
        raise NotImplementedError('Subclass to implement definition constraint')

    def evaluate_assembly(self, assembly : Assembly) -> bool:
        """Evaluate the constraint based on the complete assembled sequence."""
        raise NotImplementedError('Subclass to implement definition assembly')
