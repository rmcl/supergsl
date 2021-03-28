"""Implement the symbol table of the SuperGSL Compiler."""
from re import Pattern
from typing import Optional, Tuple, List, Dict
from supergsl.core.exception import (
    ProviderNotFoundError,
    NotFoundError,
    SymbolNotFoundError
)
from supergsl.core.provider import SuperGSLProvider
from supergsl.core.types import SuperGSLType


class SymbolTable:

    def __init__(self, name : str, parent : Optional['ScopedSymbolTable']):
        self.name = name
        self._parent = parent
        self._symbols : Dict[str, SuperGSLType] = {}
        self._nested_scopes : Dict[str, SymbolTable] = {}

    def parent_scope(self) -> Optional['ScopedSymbolTable']:
        """Return the SymbolTable for the parent scope."""
        return self._parent

    def nested_scope(self, name : str) -> 'SymbolTable':
        try:
            return self._nested_scopes[name]
        except KeyError:
            self._nested_scopes[name] = SymbolTable(name, self)
            return self._nested_scopes[name]

    def lookup(self, identifier: str) -> SuperGSLType:
        try:
            return self._symbols[identifier]
        except KeyError:
            raise SymbolNotFoundError('Symbol "{}" does not exist.'.format(identifier))

    def insert(self, identifier : str, value : SuperGSLType) -> None:
        """Set a value to identifier in this scope."""
        self._symbols[identifier] = value

    def display(self, depth=0):
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
