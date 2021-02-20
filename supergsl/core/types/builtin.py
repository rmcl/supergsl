from typing import List


class SuperGSLType(object):
    """Base class defining types available in SuperGSL."""
    pass


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
