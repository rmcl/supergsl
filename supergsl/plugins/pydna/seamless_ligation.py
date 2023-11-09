from Bio.SeqUtils import MeltingTemp as _mt
from pydna.design import assembly_fragments
from pydna.design import primer_design
from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord

from supergsl.core.assembly import AssemblerBase

from .primers import ExtractionPrimerBuilder


class SeamlessLigationAssembler(AssemblerBase):
    """Create assemblies utilizing the primer algorithm implemented by fGSL.

    Two examples of this approach are, Ligase Chain Reaction (LCR), or Seamless
    Ligation Cloning Extract (SLiCE). Check out the j5 manual for a description
    of related approaches: `"The SLIC, Gibson, CPEC, and SLiCE assembly methods"
    <https://j5.jbei.org/j5manual/pages/22.html>`_.

    .. code-block:: gsl

            from builtin import seamless_ligation

            seamless_ligation {
                HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] ; dHO
            }

    """

    import_path = 'seamless'
    name = 'seamless-ligation'

    def assemble(self, assemblies):
        primer_builder = ExtractionPrimerBuilder()

        for assembly in assemblies:
            part_records = []
            part_amplicons = []
            for part_node in assembly.parts:
                part = part_node.part

                part_amplicon, part_seq_record = primer_builder.build_primers_for_part(part)

                part_records.append(part_seq_record)
                part_amplicons.append(part_amplicon)

            fragments = assembly_fragments(part_amplicons)
            for idx in range(len(fragments)):
                fragments[idx].locus = part_records[idx].name

            #print(fragments)

            assemblyobj = Assembly(
                fragments,
                limit=20)

            linear_contigs = assemblyobj.assemble_linear()
            if len(linear_contigs) != 1:
                raise Exception(
                    '%s resulted in %d contigs. We were hoping for just one.' % (
                        len(linear_contigs)
                    ))

            assembly.contig = linear_contigs[0]
            #print(assembly.contig.figure())

        return assemblies
