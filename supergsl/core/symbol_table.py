"""Implement the symbol table of the SuperGSL Compiler."""
from typing import Optional, List, Dict, Any
from supergsl.core.exception import SymbolNotFoundError
from supergsl.types import SuperGSLType

# pylint: disable=E1136

class SymbolTable:
    """Store symbols in nested scopes."""

    def __init__(self, name : str, parent : Optional['SymbolTable']):
        self.name = name
        self._parent = parent

        # TODO: We want to restrict the type of the values of the symbol table.
        # Ideally, this would be Union[SuperGSLType, SuperGSLPlugin]. Need to resolve
        # a circular dependency first.
        self._symbols : Dict[str, Any] = {}
        self._nested_scopes : Dict[str, SymbolTable] = {}

    def parent_scope(self) -> Optional['SymbolTable']:
        """Return the SymbolTable for the parent scope."""
        return self._parent

    def enter_nested_scope(self, name : str) -> 'SymbolTable':
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
        except KeyError as error:
            raise SymbolNotFoundError(
                'Symbol "{}" does not exist.'.format(identifier)) from error

    def insert(self, identifier : str, value : Any) -> None:
        """Set a symbol to a value in the current scope."""
        self._symbols[identifier] = value

    def __iter__(self):
        """Make the table iterable returning tuples of (key, value)."""
        return iter(self._symbols.items())

    def nested_scopes(self):
        """Return a list of nested scopes"""
        return iter(self._nested_scopes.items())
