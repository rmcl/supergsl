"""Support SuperGSL's symbol provider mechanism."""
from supergsl.core.types import SuperGSLType
from supergsl.core.symbol_table import SymbolTable

class SuperGSLProvider(object):
    """Base class to define objects that can provide symbols."""

    def resolve_import(
        self,
        symbol_table : SymbolTable,
        identifier : str,
        alias : str
    ) -> None:
        """Import a identifier and register it in the symbol table."""
        raise NotImplementedError('Subclass to implement')
