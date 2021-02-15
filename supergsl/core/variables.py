import re
from re import Match, Pattern
from supergsl.core.plugin import SuperGSLPlugin

VARIABLE_MODULE_PATH = '_SuperGSLVariable'

class ResolveSymbolMixin(object):
    def resolve_import(self, identifier : str, alias : str) -> Pattern:
        """Resolve the import of a part from this provider.

        Return a tuple with:
            * A regular expression to match symbols against
            * A callback method that given the actual identifier will return the `Part`.
        """
        return re.compile(alias or identifier)

    def get_symbol(self, match : Match):
        raise NotImplementedError('Subclass to implement.')

class SuperGSLVariable(object):
    def __init__(self, identifier):
        self.identifier = identifier

        # self.type = ?????

class VariableProviderPlugin(ResolveSymbolMixin, SuperGSLPlugin):
    name = 'variable_provider'

    def register(self, symbol_table, compiler_settings):
        symbol_table.register(VARIABLE_MODULE_PATH, self)

    def get_symbol(self, match : Match):
        print(match)

        # Todo(rmcl): What do we need to return here?
        return SuperGSLVariable(match.string)
