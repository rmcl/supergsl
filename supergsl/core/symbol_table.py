"""Implement the symbol table of the SuperGSL Compiler."""
from typing import Optional, List, Dict, Any
from supergsl.core.exception import SymbolNotFoundError
from supergsl.core.types import SuperGSLType

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
        except KeyError as error:
            raise SymbolNotFoundError(
                'Symbol "{}" does not exist.'.format(identifier)) from error

    def insert(self, identifier : str, value : Any) -> None:
        """Set a symbol to a value in the current scope."""
        self._symbols[identifier] = value
