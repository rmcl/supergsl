from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.ast import SymbolRepository


class AttachSymbolRepositoryPass(BreadthFirstNodeFilteredPass):
    def __init__(self):
        self.symbol_registry = SymbolRepository()

    def visit(self, node):
        node.symbol_registry = self.symbol_registry
