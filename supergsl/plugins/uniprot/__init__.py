
from supergsl.core.plugin import SuperGSLPlugin

from .protein_provider import UniprotProvider


class UniprotPlugin(SuperGSLPlugin):
    """Plugin stub to register MixedDna SuperGSL utilities."""

    def register(self, compiler_settings):
        """Register built in assemblers."""

        self.register_available_provider(
            'uniprot',
            UniprotProvider)
