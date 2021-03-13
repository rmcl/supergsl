"""Support SuperGSL's symbol provider mechanism."""
from re import Pattern, Match
from supergsl.core.types import SuperGSLType

class SuperGSLProvider(object):
    """Base class to define objects that can provide symbols."""

    def resolve_import(self, identifier : str, alias : str) -> Pattern:
        """Resolve an identifier to be imported from this provider.

        Return a tuple with:
            * A regular expression to match symbols against
            * A callback method that given the actual identifier will return the `Part`.
        """
        raise NotImplementedError('Subclass to implement')

    def get_symbol(self, identifier_match : Match) -> SuperGSLType:
        raise NotImplementedError('Subclass to implement')
