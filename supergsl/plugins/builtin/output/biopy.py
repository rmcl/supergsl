"""Export SuperGSL Assemblies as Genbank files."""
from typing import List, Tuple, Type, TextIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from supergsl.core.types.part import Part
from supergsl.core.types.position import Slice
from supergsl.core.sequence import SequenceEntry, SequenceAnnotation
from supergsl.core.assembly import AssemblyResultOutputFunction
from supergsl.core.types.assembly import AssemblyResultSet, Assembly
from supergsl.core.constants import STRAND_WATSON, STRAND_CRICK
from supergsl.core.types.role import convert_role_to_biopython_type



def convert_to_biopy_strand(location : Slice) -> int:
    """Convert strand on `Slice` to Biopython strand."""
    if location.strand == STRAND_WATSON:
        return 1
    elif location.strand == STRAND_CRICK:
        return -1

    raise Exception('Strand not specified.')

def get_biopy_feature_type(annotation : SequenceAnnotation) -> str:
    """Try to determine the feature type for this annotation."""
    if annotation.roles:
        for role in annotation.roles:
            try:
                return convert_role_to_biopython_type(role)
            except Exception:
                pass

    if 'type' in annotation.payload:
        return annotation.payload['type']

    return 'part'

def get_feature_id_and_qualifiers(
    annotation : SequenceAnnotation,
    default_id : str
) -> Tuple[str, dict]:
    """Return the id and qualifiers for the `SeqFeature` generated from the annotation."""
    annotation_qualifiers = annotation.payload.copy()
    if annotation.roles:
        annotation_qualifiers['roles'] = [
            role.uri
            for role in annotation.roles
        ]

    possible_id_fields = ['label', 'id']
    feature_id = default_id
    for field in possible_id_fields:
        if field in annotation_qualifiers:
            feature_id = annotation_qualifiers[field]
            break

    return feature_id, annotation_qualifiers

def create_seq_feature_for_sequence_annotation(
    annotation : SequenceAnnotation,
    default_id : str,
    parent_sequence_length : int
) -> SeqFeature:
    """Create a Biopython SeqFeature for the given SequenceAnnotation."""
    start_abs_pos = annotation.location.start.build_absolute_position(
        parent_sequence_length)
    end_abs_pos = annotation.location.end.build_absolute_position(
        parent_sequence_length)

    feature_id, qualifiers = get_feature_id_and_qualifiers(annotation, default_id)
    if annotation.location.strand == STRAND_CRICK:
        tmp_end_abs_pos = start_abs_pos.get_complement_strand_position()
        start_abs_pos = end_abs_pos.get_complement_strand_position()
        end_abs_pos = tmp_end_abs_pos

    feature = SeqFeature(
        id=feature_id,
        qualifiers=qualifiers,
        location=FeatureLocation(
            start=start_abs_pos.index,
            end=end_abs_pos.index
        ),
        strand=convert_to_biopy_strand(annotation.location),
        type=get_biopy_feature_type(annotation))

    return feature

def create_seq_record_from_sequence_entry(
    sequence_entry : SequenceEntry,
    record_name : str
) -> SeqRecord:
    """Create a Biopython SeqRecord for a SequenceEntry"""
    record = SeqRecord(
        sequence_entry.sequence,
        id=str(sequence_entry.id),
        name=record_name,
        description='')

    for annotation_idx, annotation in enumerate(sequence_entry.annotations()):
        feature = create_seq_feature_for_sequence_annotation(
            annotation,
            str(annotation_idx),
            sequence_entry.sequence_length)

        record.features.append(feature)

    return record

def build_seq_record_for_assembly(assembly : Assembly) -> SeqRecord:
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

    sequence_entry = assembly.part.sequence_entry
    for annotation_idx, annotation in enumerate(sequence_entry.annotations()):
        feature = create_seq_feature_for_sequence_annotation(
            annotation,
            str(annotation_idx),
            sequence_entry.sequence_length)

        record.features.append(feature)

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



class SeqRecordAssemblyOutput(AssemblyResultOutputFunction):
    """Generate file containing an annotated assembly.

    Output format defaults to genbank, but can be anything supported by BioPython.
    """

    def output(self, assemblies : AssemblyResultSet, file_handle : TextIO, **kwargs):
        records : List[SeqRecord] = []

        file_format = kwargs.get('file_format', 'genbank')
        for assembly in assemblies:
            records.append(
                build_seq_record_for_assembly(assembly))

        SeqIO.write(records, file_handle, file_format)
