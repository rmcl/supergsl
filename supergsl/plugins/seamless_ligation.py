from supergsl.core.assembly import AssemblerBase

from pydna.design import assembly_fragments
from pydna.design import primer_design
from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord

class SeamlessLigationAssembler(AssemblerBase):
    """Create assemblies utilizing the primer algorithm implemented by fGSL."""

    @property
    def target_Tm(self):
        try:
            return self.options['target_Tm']
        except KeyError:
            raise Exception('`target_Tm` must be specified in Assembler options.')

    def assemble(self, assemblies):

        for assembly in assemblies:

            part_records = []
            part_amplicons = []
            for part_node in assembly.parts:
                part = part_node.part
                part_seq_record = part.get_sequence()

                part_has_primers = (
                    True if part.forward_primer and part.reverse_primer
                    else False
                )

                amplicon = primer_design(
                    Dseqrecord(part_seq_record),
                    fp=part.forward_primer,
                    rp=part.reverse_primer,
                    target_tm=self.target_Tm,
                    limit=13)

                if not part_has_primers:
                    part.set_primers(
                        amplicon.forward_primer,
                        amplicon.reverse_primer)

                part_records.append(part_seq_record)
                part_amplicons.append(amplicon)

            fragments = assembly_fragments(part_amplicons)
            for idx in range(len(fragments)):
                fragments[idx].locus = part_records[idx].name

            print(fragments)

            assemblyobj = Assembly(
                fragments,
                limit=20)

            lin_asms = assemblyobj.assemble_linear()
            linear = lin_asms[0]
            #print(len(lin_asms), lin_asms)
            #print(assemblyobj)
            print(linear.figure())

            #import pdb
            #pdb.set_trace()
