from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunctionDeclaration
from .biobrick import BioBrick3AAssembler, BioBrickPartProvider

class BioBrickPlugin(SuperGSLPlugin):

    def register(self, compiler_settings):
        """Register functions provide by the igem plugin."""

        self.register_function(
            'biobrick',
            'assemble-3a',
            SuperGSLFunctionDeclaration(BioBrick3AAssembler, compiler_settings))

        self.register_provider(
            'biobrick',
            BioBrickPartProvider('biobrick', compiler_settings))
