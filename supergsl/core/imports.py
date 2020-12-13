from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.exception import (
    ProviderNotFoundException,
    FunctionNotFoundException,
    GSLImportError
)


class ResolveImportsPass(BreadthFirstNodeFilteredPass):

    def get_node_handlers(self):
        return {
            'ProgramImport': self.visit_import_node,
            'Part': self.visit_part_node,
            'FunctionInvocation': self.visit_function_invoke_node,
        }

    def before_pass(self, ast):
        return ast

    def resolve_import(self, module_path, identifier):
        part_symbol_table = self.symbol_registry.get_table('parts')

        try:
            part_symbol_table.resolve_part(module_path, identifier)
            return
        except ProviderNotFoundException:
            pass

        function_symbol_table = self.symbol_registry.get_table('functions')
        try:
            function_symbol_table.resolve_function(module_path, identifier)
        except FunctionNotFoundException:
            raise GSLImportError('"%s" did not contain "%s".' % (
                module_path,
                identifier
            ))


    def visit_import_node(self, node):
        for program_import in node.imports:
            module_path = '.'.join(node.module)
            self.resolve_import(module_path, program_import.identifier)

        return node

    def visit_function_invoke_node(self, node):
        function_symbol_table = self.symbol_registry.get_table('functions')
        node.function = function_symbol_table.get_function(node.identifier)
        return node

    def visit_part_node(self, node):
        part_symbol_table = self.symbol_registry.get_table('parts')
        node.part = part_symbol_table.get_part(node.identifier)
        return node
