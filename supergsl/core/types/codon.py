from .base import SuperGSLType


class CodonFrequencyTable(SuperGSLType):
    """Store a table of codon frequencies.

    https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables/tree/master/codon_usage_data/tables
    """
    def __init__(self, name : str, table : dict):
        self.name = name
        self.table = table

    def __repr__(self):
        return self.name

    def print(self):
        """Display the table of amino acid, codon, and frequency."""
        result = [
            'Codon Frequencies for: %s' % self.name,
            'AmAcid\tCodon\tFrequency'
        ]
        for amino_acid, codons in self.table.items():
            for codon, frequency in codons.items():
                result.append('%s\t%s\t%06f' % (
                    amino_acid, codon, frequency
                ))
            result.append('')

        return '\n'.join(result)
