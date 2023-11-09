"""Define a function for instantiating AssemblyDeclarations."""
from typing import List
from supergsl.core.types.assembly import AssemblyDeclaration
from supergsl.core.function import SuperGSLFunction, SuperGSLFunctionDeclaration


class AssemblyDeclarationFunction(SuperGSLFunction):
    """Define a function to allow for assembly declarations.

    ```
    let x = declare {
        part1 ; part2 ; part3
    }
    ```
    """
    return_type = List[AssemblyDeclaration]

    def execute(self, params):
        return params['children']
