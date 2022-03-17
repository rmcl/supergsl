"""Define the Assembler base class."""
import sys
import contextlib
from typing import List, TextIO

from supergsl.core.function import SuperGSLFunction
from supergsl.core.types.assembly import (
    AssemblyDeclaration,
    AssemblyResultSet,
    Assembly,
    AssemblyLevel
)
from supergsl.core.types.constraint import Constraint


class AssemblerBase(SuperGSLFunction):
    """Base class for functions implementing Assemblers."""

    return_type = AssemblyResultSet

    def execute(self, params):
        assembly_requests : List[AssemblyDeclaration] = params['children']
        contraints : List[Constraint] = params['constraints']

        result_set = AssemblyResultSet([])
        for assembly_idx, assembly_request in enumerate(assembly_requests):

            designs = assembly_request.get_designs()
            for design_idx, design_description in enumerate(designs):

                if not self.evaluate_definition_contraints(design_description, constraints):
                    continue

                assembly = self.assemble(design_idx, design_description)

                if not self.evaluate_sequence_constraints(assembly, constraints):
                    continue

                result_set.add_assembly(assembly)

        return result_set

    def evaluate_definition_contraints(
        self,
        design_description : List[AssemblyLevel],
        constraints : List[Constraint]
    ) -> bool:
        """Evaluate constraints based on just the definition of a design.

        This is executed *before* the assembler runs and allows the compiler to skip
        designs that can be removed from consideration based only on the definition.
        """
        for constraint in constraints:
            if not constraint.is_definition_constraint:
                return True

            return constraint.evaluate_definition(design_description)

    def evaluate_assembly_constraints(
        self,
        assembly : Assembly,
        constraints : List[Constraint]
    ) -> bool:
        """Evaluate constraints based on the complete assembly of a design.

        This is executed *after* the assembler runs and should give the constraint
        full access to the complete sequence of the assembly.
        """
        for constraint in constraints:
            if not constraint.is_assembly_constraint:
                return True

            return constraint.evaluate_assembly(assembly)

    def assemble(self, assembly : List[AssemblyLevel]) -> Assembly:
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
