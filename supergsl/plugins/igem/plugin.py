from supergsl.core.plugin import SuperGSLPlugin

from .oligos import SyntheticOligoAssembler
from .biobrick import BioBrick3AAssembler, BioBrickPartProvider

class BioBrickPlugin(SuperGSLPlugin):

    def register(self, symbol_table, compiler_settings):
        """Register functions provide by the igem plugin."""

        symbol_table.register('synthesis', SyntheticOligoAssembler(compiler_settings))
        symbol_table.register('biobrick', BioBrick3AAssembler(compiler_settings))
        symbol_table.register('biobrick', BioBrickPartProvider('biobrick', compiler_settings))
