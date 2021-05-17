from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunctionDeclaration

from .provider import CodonFrequencyTableProvider

class CodonFrequencyTablePlugin(SuperGSLPlugin):
    """Plugin stub to register Codon Frequecy Table loader."""

    def register(self, compiler_settings : dict):
        """Register frequency table provider."""
        self.register_provider(
            'codon_frequency',
            CodonFrequencyTableProvider(compiler_settings))
