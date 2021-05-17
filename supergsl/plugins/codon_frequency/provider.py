from typing import Dict, Tuple, List

from supergsl.core.types.builtin import CodonFrequencyTable
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.provider import SuperGSLProvider

from python_codon_tables import (
    available_codon_tables_names,
    get_codons_table
)


class CodonFrequencyTableProvider(SuperGSLProvider):
    """Retrieve codon usage tables for popular organisms.

    All the tables are from kazusa.or.jp and the original paper:

    Codon usage tabulated from the international DNA sequence databases:
    status for the year 2000.
    Nakamura, Y., Gojobori, T. and Ikemura, T. (2000) Nucl. Acids Res. 28, 292.

    To configure this part provider add the following to `supergsl-config.json`:
    ```
    {
        "name": "kazusa_codon_tables",
        "provider_class": "supergsl.plugins.codon_frequency.CodonFrequencyTableProvider"
    }

    """

    def __init__(self, compiler_settings : dict):
        self.settings = compiler_settings

    def list_tables(self):
        return available_codon_tables_names

    def get_table(self, identifier : str) -> CodonFrequencyTable:
        """Retrieve a codon frequency table by organism identifier.

        Arguments:
            identifier  A identifier to select a codon table
        Return: `CodonFrequencyTable`
        """
        table = get_codons_table(identifier)

        # TODO: Figure out how we convert this table to a gsl type.
        return CodonFrequencyTable()

    def resolve_import(
        self,
        symbol_table : SymbolTable,
        identifier : str,
        alias : str
    ) -> None:
        """Import a table and register it in the symbol table."""
        table_identifier = alias or identifier
        symbol_table.insert(table_identifier, self.get_table(identifier))
