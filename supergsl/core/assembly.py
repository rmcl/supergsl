from typing import List
from supergsl.core.function import SuperGSLFunction
from supergsl.core.types.assembly import AssemblyDeclaration, AssemblyList


class AssemblerBase(SuperGSLFunction):
    """Base class for functions implementing Assemblers."""
    def get_arguments(self):
        return []

    def get_return_type(self):
        return AssemblyList

    def execute(self, params):
        return self.assemble(params['children'])

    def assemble(self, assembly_requests : List[AssemblyDeclaration]) -> AssemblyList:
        """Iterate over `Part` and generate an Assembly object."""
        raise NotImplementedError('Not implemented. Subclass to implement.')
