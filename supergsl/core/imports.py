from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.exception import (
    ProviderNotFoundException,
    FunctionNotFoundException,
    GSLImportError
)

from supergsl.core.variables import VARIABLE_MODULE_PATH


class ResolveImportsPass(BreadthFirstNodeFilteredPass):

    def get_node_handlers(self):
        return {
            'Import': self.visit_import_node,
            'Part': self.visit_part_node,
            'FunctionInvocation': self.visit_function_invoke_node,
            'VariableDeclaration': self.visit_variable_declaration_node,
        }

    def before_pass(self, ast):
        return ast

    def visit_import_node(self, node):
        for program_import in node.imports:
            module_path = '.'.join(node.module)
            node.symbol = self.symbol_table.resolve_symbol(
                module_path,
                program_import.identifier)

        return node

    def visit_variable_declaration_node(self, node):
        node.variable = self.symbol_table.resolve_symbol(
            VARIABLE_MODULE_PATH,
            node.identifier)
        return node

    def visit_function_invoke_node(self, node):
        node.function = self.symbol_table.get_symbol(node.identifier)

        ## TODO: DO TYPE CHECKING TO ASSERT THIS SYMBOL IS A FUNCTION

        return node

    def visit_part_node(self, node):
        node.part = self.symbol_table.get_symbol(node.identifier)
        return node
