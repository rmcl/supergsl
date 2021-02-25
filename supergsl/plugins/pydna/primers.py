"""Generate Extraction Primers using the excellent pydna wrapper around BioPython."""

from typing import Callable, Tuple
from Bio.SeqUtils import MeltingTemp as _mt
from pydna.design import assembly_fragments
from pydna.design import primer_design
from pydna.dseqrecord import Dseqrecord

from supergsl.core.types import PrimerPair


class ExtractionPrimerBuilder(object):
    """Use pydna to generate extraction primers for parts."""

    def __init__(self, options = None):
        self.options = options or {}

    def build_primers_for_part(self, part):
        """Build extraction primers and assign them into the part."""
        part_seq_record = part.get_sequence()

        forward_primer = part.extraction_primers.forward_primer if part.has_primers else None
        reverse_primer = part.extraction_primers.reverse_primer if part.has_primers else None

        amplicon = primer_design(
            Dseqrecord(part_seq_record),
            fp=forward_primer,
            rp=reverse_primer,
            tm_func=self._perform_tm_func,
            limit=13)

        if not part.has_primers:
            part.set_extraction_primers(
                PrimerPair.from_sequences(
                    amplicon.forward_primer.seq,
                    amplicon.reverse_primer.seq
                ))

        return amplicon

    def _perform_tm_func(self, seq):
        """Call the tm_func on a sequence with the config specified arguments."""
        tm_func, tm_kwargs = self._tm_func_parameters()
        return tm_func(seq, **tm_kwargs)

    def _tm_func_parameters(self) -> Tuple[Callable, dict]:
        """Override the BioPython tm calculation method.

        These defaults have largely been pulled from pyDNA, however here we
        choose to use DNA_NN3 because that is the nucleotide melting table
        used by fGSL [2].

        A note from pyDNA source [1]: See the documentation for Bio.SeqUtils.MeltingTemp for more details
        The 10X Taq Buffer with (NH4)2SO4 is commercialized by companies like
            * ThermoFisher, although we make it ourselves
            * 10X Buffer Composition
            * 750 mM Tris-HCl (pH 8.8 at 25°C),
            * 200 mM (NH4)2SO4,
            * 0.1% (v/v) Tween 20.

        References:
        1. https://github.com/BjornFJohansson/pydna/blob/master/src/pydna/tm.py#L17
        2. https://github.com/Amyris/AmyrisBio/blob/master/src/AmyrisBio/primercore2.fs#L173
        """

        if 'tm_func' not in self.options:
            # Used by Primer3Plus to calculate the product Tm.
            tm_func = _mt.Tm_NN
        else:
            raise Exception('Need to implement swapping out default Tm_NN function.')

        tm_kwargs = {
            'check': True,
            'strict': True,
            'c_seq': None,
            'shift': 0,

            # DNA_NN4: values from SantaLucia & Hicks (2004) is the default used by pydna
            # however DNA_NN3 is effectively the value used by the fGSL primer designer.
            'nn_table': _mt.DNA_NN3,
            'tmm_table': None,
            'imm_table': None,
            'de_table': None,
            'dnac1': 500 / 2,  # Assume 500 µM of each primer in the PCR mix
            'dnac2': 500 / 2,  # This is what MELTING and Primer3Plus do
            'selfcomp': False,
            'Na': 40,
            'K': 0,
            'Tris': 75.0,  # We use the 10X Taq Buffer with (NH4)2SO4 (above)
            'Mg': 1.5,  # 1.5 mM Mg2+ is often seen in modern protocols
            'dNTPs': 0.8,  # Assume 200 µM of each dNTP
            'saltcorr': 7,  # Tm = 81.5 + 0.41(%GC) - 600/N + 16.6 x log[Na+]
        }

        for key in tm_kwargs.keys():
            if key in self.options:
                tm_kwargs[key] = self.options[key]

        return (tm_func, tm_kwargs)
