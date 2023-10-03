from Bio.Seq import Seq

# From the docs, "BioBrick RFC[10] it must not contain the following
# restriction sites, as these are unique to the prefix and suffix"
# Taken from http://parts.igem.org/Help:Standards/Assembly/RFC10
BIOBRICK_ILLEGAL_RESTRICTION_SITES = ['EcoRI', 'XbaI', 'SpeI', 'PstI', 'NotI']

BIOBRICK_PREFIX_CODING_SEQUENCE = Seq('GAATTCGCGGCCGCTTCTAG')
BIOBRICK_PREFIX_SEQUENCE = Seq('GAATTCGCGGCCGCTTCTAGAG')
BIOBRICK_SUFFIX_SEQUENCE = Seq('TACTAGTAGCGGCCGCTGCAG')
