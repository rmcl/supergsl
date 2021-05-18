"""Retrieve codon usage tables for popular organisms.

All the tables are from kazusa.or.jp and the original paper:

Codon usage tabulated from the international DNA sequence databases:
status for the year 2000.
Nakamura, Y., Gojobori, T. and Ikemura, T. (2000) Nucl. Acids Res. 28, 292.

To configure this plugin add the following path to plugins in `supergsl-config.json`:
```
{
    ...
    "plugins": [
        ...
        "supergsl.plugins.codon_frequency"
    ]
}
```
"""
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
