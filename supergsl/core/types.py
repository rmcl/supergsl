class SuperGSLType(object):
    """Base class defining types available in SuperGSL."""
    pass


class SuperGSLEnum(SuperGSLType):
    """Define a list of choices."""
    options = []


class NucleotideSequence(SuperGSLType):
    """A type representing a nucleotide sequence."""

    def get_sequence(self):
        """Return the nucleotide sequence as a `Bio.Seq`."""
        raise NotImplementedError('Subclass to implement.')


class AminoAcidSequence(SuperGSLType):
    """A type representing arbitrary an amino acid sequence."""

    def get_sequence(self):
        """Return the amino acid sequence as a `Bio.Seq`."""
        raise NotImplementedError('Subclass to implement.')


class CodonFrequencyTable(SuperGSLType):
    """
    https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables/tree/master/codon_usage_data/tables
    """


class CodonTranslationTable(SuperGSLType):
    """Encode the mapping of DNA/RNA codons to protein.

    This is a thin wrapper around biopython's `Bio.Data.CodonTable`
    https://biopython.org/docs/1.75/api/Bio.Data.CodonTable.html
    """
