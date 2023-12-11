
from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.provider import ProviderConfig
from .protein_provider import UniprotProvider


class UniprotPlugin(SuperGSLPlugin):
    """Plugin stub to register MixedDna SuperGSL utilities."""

    def register(self, compiler_settings):
        """Register built in assemblers."""

        self.register_available_provider(
            'uniprot',
            UniprotProvider)

        # Todo: We gotta improve the interface to creating providers config.
        sequence_store = self.symbol_table.lookup('sequences')
        config = ProviderConfig(
            sequence_store, compiler_settings)

        self.register_provider(
            'uniprot',
            UniprotProvider('uniprot', config))
