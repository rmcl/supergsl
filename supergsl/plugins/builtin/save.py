"""Simple function to save SuperGSLTypes to a provider."""
from supergsl.core.types import SuperGSLType
from supergsl.core.provider import SuperGSLProvider
from supergsl.core.function import SuperGSLFunction
from supergsl.core.types.builtin import Collection


class SaveFunction(SuperGSLFunction):
    """A simple function to save a part.

    This is primarily useful in the REPL shell.
    """
    arguments = [
        ('provider', SuperGSLProvider),
        ('item', SuperGSLType)
    ]

    def execute(self, params : dict):
        """Save a part to the provider."""

        if isinstance(params['item'], Collection):
            for item in params['item']:
                params['provider'].save(item)
        else:
            params['provider'].save(params['item'])
