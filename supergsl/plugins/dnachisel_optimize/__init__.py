from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunctionDeclaration

from .function import DNAChiselOptimizeFunction

class DNAChiselPlugin(SuperGSLPlugin):
    """Plugin stub to help register basic Assemblers."""

    def register(self, compiler_settings : dict):
        """Register built in assemblers."""
        self.register_function(
            'dnachisel',
            'codon_optimize',
            SuperGSLFunctionDeclaration(
                DNAChiselOptimizeFunction,
                compiler_settings))
