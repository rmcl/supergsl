from Bio.Seq import Seq
from supergsl.core.types.assembly import Assembly
from supergsl.core.assembly import AssemblerBase

DEFAULT_MAX_OLIGO_LEN = 200
DEFAULT_MIN_OVERLAP = 20
DEFAULT_NUM_OLIGOS = 2


class SyntheticOligoAssembler(AssemblerBase):
    """Create parts that to constructed by ordering a set of short synthetic oligos.

    These oligos are typically used as primers, but can be used for short parts
    without a template.
    """

    import_path = 'basic'
    name = 'synthetic-oligos'

    def __init__(self, config_options):
        self.max_oligo_len = config_options.get('max_oligo_len', DEFAULT_MAX_OLIGO_LEN)
        self.min_overlap = config_options.get('max_oligo_len', DEFAULT_MIN_OVERLAP)
        self.max_num_oligos = config_options.get('max_num_oligos', DEFAULT_NUM_OLIGOS)


    def assemble(self, assemblies):
        """Assemble for synthesis using synthetic oligos."""
        oligos = []
        for assembly_node in assemblies.definitions:
            part_sequences = []
            for part_node in assembly_node.parts:
                part = part_node.part

                print('PART', part)
                part_sequences.append(part.get_sequence().seq)

            # Todo: Find a better way to concatenate Bio.Seq objects.
            assembly_sequence = Seq(''.join(str(part_sequences)))

            oligo_idx = 0
            seq_pos = 0
            sequence_length = len(assembly_sequence)
            while True:
                if oligo_idx == self.max_num_oligos:
                    raise Exception((
                        'Assembly is to large to be synthesized. Construct length: {} '
                        'Max Oligo Length: {} Max Oligos: {} ').format(
                            len(assembly_sequence),
                            self.max_oligo_len,
                            self.max_num_oligos
                        ))

                remaining_seq_len = sequence_length - seq_pos
                if remaining_seq_len < self.max_oligo_len:
                    seq_pos = min(0, seq_pos - (self.max_oligo_len - remaining_seq_len))

                oligos.append(assembly_sequence[seq_pos:self.max_oligo_len])
                seq_pos += self.max_oligo_len
                if seq_pos >= sequence_length:
                    break

                seq_pos -= self.min_overlap
                oligo_idx += 1

            print(oligos)
            new_assembly = Assembly(assembly_sequence, oligos)
