"""Define the Assembler base class."""
import sys
import contextlib
from typing import List, TextIO

from supergsl.core.function import SuperGSLFunction
from supergsl.core.types.assembly import (
    AssemblyDeclaration,
    AssemblyResultSet,
    AssemblyLevel,
    Assembly
)


class AssemblerBase(SuperGSLFunction):
    """Base class for functions implementing Assemblers."""

    return_type = AssemblyResultSet

    def execute(self, params : dict):
        assembly_requests : List[AssemblyDeclaration] = params['children']

        return self.assemble(assembly_requests)

    def assemble(
        self,
        assembly_requests : List[AssemblyDeclaration]
    ) -> AssemblyResultSet:
        """Assemble a set of AssemblyDeclarations to create an AssemblyResultSet."""
        result_set = AssemblyResultSet([])
        for assembly_idx, assembly_request in enumerate(assembly_requests):
            assembly_label = assembly_request.label or f'{assembly_idx:03}'

            designs = assembly_request.get_designs()
            for design_idx, design_description in enumerate(designs):
                design_label = f'{design_idx:03}'
                assembly = self.assemble_design(assembly_label, design_label, design_description)
                result_set.add_assembly(assembly)

        return result_set


    def assemble_design(
        self,
        assembly_label : str,
        design_label : str,
        assembly_request : List[AssemblyLevel]
    ) -> Assembly:
        """Iterate over a list of `Part` and generate an Assembly object."""
        raise NotImplementedError('Not implemented. Subclass to implement.')


class AssemblyResultOutputFunction(SuperGSLFunction):
    """Base SuperGSLFunction for creating exporters for AssemblyResultSets.

    When subclassed creates a SuperGSL Function with the following parameters:
    ```output_sbol(assemblies, filename)``` and the side effect of creating a
    file (filename) with outputted AssemblyResultSet.

    Good examples:
        supergsl/plugins/builtin/output/json_output.py
        supergsl/plugins/builtin/output/sbol_output.py

    """

    def output(self, assemblies : AssemblyResultSet, file_handle : TextIO):
        """Output an `AssemblyResultSet` to the supplied file_handle."""
        raise NotImplementedError('Subclass to implement.')

    def get_arguments(self):
        return [
            ('assemblies', AssemblyResultSet),
            ('filename', str),
        ]

    def get_return_type(self):
        return type(None)

    @contextlib.contextmanager
    def open_output_fp(self, filename = None):
        """Open the output file. Use contextmanager to do smart stuff around stdout."""
        if filename and filename != '-':
            file_handle = open(filename, 'w+')
        else:
            file_handle = sys.stdout

        try:
            yield file_handle
        finally:
            if file_handle is not sys.stdout:
                file_handle.close()

    def execute(self, params : dict):
        """Execute the Output Function."""
        with self.open_output_fp(params['filename']) as output_fp:
            self.output(params['assemblies'], output_fp)
