from supergsl.core.plugin import SuperGSLPlugin

from .part_library import MixedPartLibraryProvider


class MixedDnaPlugin(SuperGSLPlugin):
    """Plugin stub to register MixedDna SuperGSL utilities."""

    def register(self, compiler_settings):
        """Register built in assemblers."""

        self.register_available_provider(
            'mixed_part_library',
            MixedPartLibraryProvider)
