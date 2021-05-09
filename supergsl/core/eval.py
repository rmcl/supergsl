"""Define the mechanism of SuperGSLFunction and an AST pass to invoke those functions."""
from typing import List, Dict, Set

from supergsl.core.types import SuperGSLType
from supergsl.core.backend import DepthFirstNodeFilteredPass
from supergsl.core.types.builtin import Collection
from supergsl.utils import display_symbol_table
from supergsl.core.ast import Program as AstProgram

from supergsl.core.exception import FunctionInvokeError

#pylint: disable=E1136


class EvaluatePass(DepthFirstNodeFilteredPass):
    """Traverse the AST to build a depenency graph. Then execute all nodes in dependency order."""

    def get_node_handlers(self):
        return {
            'Program': self.visit_program_node,
        }

    def before_pass(self, ast):
        #display_symbol_table(self.symbol_table)

        self.unassigned_assembly_counter : int = 0
        self.unassigned_symbol_counter : int = 0
        return ast

    def after_pass(self, ast):
        display_symbol_table(self.symbol_table)
        return ast

    def visit_program_node(self, program_node : AstProgram):

        if not program_node.definitions:
            print('WARNING NO DEFINITIONS!!')

        symbol_table = self.get_symbol_table()
        for definition_node in program_node.definitions.definitions:
            node_type = type(definition_node).__name__

            if node_type == 'VariableDeclaration':
                definition_node.eval()
            elif node_type == 'Assembly':
                self.unassigned_symbol_counter += 1
                result = definition_node.eval()

                temp_identifer = 'assembly_%03d' % self.unassigned_symbol_counter
                symbol_table.insert(temp_identifer, result)

            elif node_type == 'FunctionInvocation':
                self.unassigned_symbol_counter += 1
                temp_identifer = 'func_%03d' % self.unassigned_symbol_counter
                result = definition_node.eval()
                symbol_table.insert(temp_identifer, result)

        return program_node
