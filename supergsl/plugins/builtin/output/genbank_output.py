"""Export SuperGSL Assemblies as Genbank files."""
from typing import List, Tuple, Type, TextIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from supergsl.core.types.part import Part
from supergsl.core.assembly import AssemblyResultOutputFunction
from supergsl.core.types.assembly import AssemblyResultSet, Assembly


class GenBankOutput(AssemblyResultOutputFunction):
    """Generate GenBank file containing an annotated assembly."""

    def output(self, assemblies : AssemblyResultSet, file_handle : TextIO):
        records : List[SeqRecord] = []

        for assembly in assemblies:
            records.append(
                self.build_seq_record_for_assembly(assembly))

        SeqIO.write(records, file_handle, 'genbank')

    def build_seq_record_for_assembly(self, assembly : Assembly) -> SeqRecord:
        """Build a `SeqRecord` entry for the given `Assembly`."""

        assembly_sequence = assembly.part.sequence
        record = SeqRecord(
            assembly_sequence,
            id='123456789', # random accession number
            name=assembly.identifier,
            description=assembly.description or '')


        sequence_entry_to_part = {}
        for reagent in assembly.reagents:
            if not isinstance(reagent, Part):
                # Exclude all reagents that are not parts
                continue

            sequence_entry_to_part[reagent.sequence_entry.id] = reagent

        for parent_link in assembly.part.sequence_entry.parent_links:
            parent_entry = parent_link.parent_entry
            target_slice = parent_link.target_slice
            start_abs_pos = target_slice.start.build_absolute_position(len(assembly_sequence))
            end_abs_pos = target_slice.end.build_absolute_position(len(assembly_sequence))

            parent_part = sequence_entry_to_part[parent_entry.id]

            feature = SeqFeature(
                id=parent_part.identifier,
                qualifiers={
                    'name': parent_part.identifier,
                    'description': parent_part.description,
                    'roles': [
                        role.uri
                        for role in parent_link.roles
                    ]
                },
                location=FeatureLocation(
                    start=start_abs_pos.index,
                    end=end_abs_pos.index
                ),
                type='part')

            record.features.append(feature)

        return record
