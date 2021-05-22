from supergsl.core.function import SuperGSLFunctionDeclaration
from supergsl.core.plugin import SuperGSLPlugin

from .cut_cas9 import GuideRNABuilderFunction


class ChopChopPlugin(SuperGSLPlugin):
    def register(self, compiler_settings : dict):
        """Register functions provide by chopchop."""
        print('REGISTER CHOP')
        self.register_function(
            'chopchop',
            'cut_cas9',
            SuperGSLFunctionDeclaration(
                GuideRNABuilderFunction, compiler_settings))
