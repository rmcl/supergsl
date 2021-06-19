from Bio import SeqIO
from Bio.Restriction import RestrictionBatch

from .exception import InvalidBioBrickError
from .constants import BIOBRICK_ILLEGAL_RESTRICTION_SITES


def check_is_valid_biobrick(sequence):
    """A valid biobrick part does not have any disallowed restriction sites except in
    the prefix and suffix.

    Return True if sequence is valid or raise `InvalidSequenceError`.
    """
    illegal_sites = RestrictionBatch(BIOBRICK_ILLEGAL_RESTRICTION_SITES)
    results = illegal_sites.search(sequence)

    for restriction_site, match_positions in results.items():

        if len(match_positions) > 0:
            raise InvalidBioBrickError(
                '{} site found in sequence.'.format(restriction_site))

    return True

def load_biobrick_constant_sequences():
    file_path = 'supergsl/plugins/igem/biobrick.fa'
    records = SeqIO.parse(file_path, "fasta")

    return {
        record.name: record
        for record in records
    }
