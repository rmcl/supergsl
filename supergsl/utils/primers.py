from typing import Dict

from Bio.Seq import Seq
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design, assembly_fragments
from pydna.assembly import Assembly

from supergsl.core.types.primer import Primer

def build_two_part_assembly_primers(
    five_prime_seq: Seq,
    three_prime_seq : Seq,
    payload_seq : Seq
) -> Dict[str, Seq]:
    """Select primers to assemble two sequences with an optional payload in-between.

    *****FIVE-PRIME SEQ*****
                            **PAYLOAD**
                                       *****THREE-PRIME SEQ*****
    terminal forward --->
                                        <<<--- terminal reverse
                        internal forward-->
                      <---internal reverse

    """

    five_prime_amplicon = primer_design(Dseqrecord(five_prime_seq))
    three_prime_amplicon = primer_design(Dseqrecord(three_prime_seq))

    fp_fragment, tp_fragment = assembly_fragments([
        five_prime_amplicon,
        Dseqrecord(payload_seq),
        three_prime_amplicon
    ])

    assembly = Assembly([fp_fragment, tp_fragment])
    contig = assembly.assemble_linear()[0]

    return {
        'terminal_forward_primer': Primer(fp_fragment.forward_primer.seq.reverse_complement()),
        'terminal_reverse_primer': Primer(tp_fragment.reverse_primer.seq),

        'internal_forward_primer': Primer(tp_fragment.forward_primer.seq.reverse_complement()),
        'internal_reverse_primer': Primer(fp_fragment.reverse_primer.seq),

        'assembly_sequence': Seq(contig.seq.watson)
    }
