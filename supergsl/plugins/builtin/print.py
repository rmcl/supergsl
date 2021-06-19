"""Simple function to print SuperGSLTypes."""
from supergsl.core.types import SuperGSLType
from supergsl.core.function import SuperGSLFunction


class SuperGSLTypeDetailFunction(SuperGSLFunction):
    """A simple function to print details of a SuperGSLTypes.

    This is primarily useful in the REPL shell.
    """
    arguments = [
        ('item', SuperGSLType)
    ]

    def execute(self, params : dict):
        print(params['item'].print())


class SuperGSLTypeHelpFunction(SuperGSLFunction):
    """A simple function to print help from SuperGSLTypes.

    This is primarily useful in the REPL shell.
    """

    arguments = [
        ('item', SuperGSLType)
    ]

    def execute(self, params : dict):
        print(params['item'].help)
