"""Simple function to print SuperGSLTypes."""
from supergsl.types import SuperGSLType
from supergsl.core.function import SuperGSLFunction


class SuperGSLTypePrintFunction(SuperGSLFunction):
    """A simple function to print SuperGSLTypes.

    This is primarily useful in the REPL shell.
    """

    def get_arguments(self):
        return [
            ('item', SuperGSLType)
        ]

    def get_return_type(self):
        return None

    def execute(self, params : dict):
        print(params['item'].print())
