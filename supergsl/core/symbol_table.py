"""Implement the symbol table of the SuperGSL Compiler."""
from typing import Optional, List, Dict
from supergsl.core.exception import SymbolNotFoundError
from supergsl.core.types import SuperGSLType

# pylint: disable=E1136

class SymbolTable:
    """Store symbols in nested scopes."""

    def __init__(self, name : str, parent : Optional['SymbolTable']):
        self.name = name
        self._parent = parent
        self._symbols : Dict[str, SuperGSLType] = {}
        self._nested_scopes : Dict[str, SymbolTable] = {}

    def parent_scope(self) -> Optional['SymbolTable']:
        """Return the SymbolTable for the parent scope."""
        return self._parent

    def nested_scope(self, name : str) -> 'SymbolTable':
        """Get or Create a nested scope within the current table."""
        try:
            return self._nested_scopes[name]
        except KeyError:
            self._nested_scopes[name] = SymbolTable(name, self)
            return self._nested_scopes[name]

    def lookup(self, identifier: str) -> SuperGSLType:
        """Lookup a symbol in the table."""
        try:
            return self._symbols[identifier]
        except KeyError:
            raise SymbolNotFoundError('Symbol "{}" does not exist.'.format(identifier))

    def insert(self, identifier : str, value : SuperGSLType) -> None:
        """Set a symbol to a value in the current scope."""
        self._symbols[identifier] = value

    def display(self, depth=0):
        """Recursively display the table and all of its nested scopes."""
        indent = ' ' * (depth*4)
        print(indent + 'Table: {} ----'.format(self.name))
        for key, val in self._symbols.items():
            print(indent + '  Key: {}, Type: {}'.format(
                key,
                type(val)
            ))

        for nested_scope in self._nested_scopes.values():
            nested_scope.display(depth+1)

        print(indent + '--- End Table: %s ----')
