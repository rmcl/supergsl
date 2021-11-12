"""Functions for serializing superGSL objects for gRPC."""
from supergsl.core.sequence import (
    SequenceEntry,
    EntryLink,
    Role
)
from supergsl.grpc.stubs.sgsl_pb2 import (
    SequenceEntry as gRPCSequenceEntry,
    SequenceLink,
    SequenceSlice,
    Role as gRPCRole
)

def role_dto(
    role : Role
) -> gRPCRole:
    """Serialize a SuperGSL Sequence Role as a gRPC Role."""
    return gRPCRole(
        uri=role.uri,
        name=role.name,
        description=role.description)

def entry_link_dto(
    entry_link : EntryLink
) -> SequenceLink:
    """Serialize a EntryLink as a gRPC SequenceLink."""
    sequence_link = SequenceLink(
        parent_entry_identifier=str(entry_link.parent_entry.id),
        child_entry_identifier=str(sequence_entry.id),
        source_slice=SequenceSlice(),
        target_slice=SequenceSlice())

    sequence_link.roles.extend([
        role_dto(role)
        for role in entry_link
    ])
    return sequence_link

def sequence_entry_dto(
    sequence_entry : SequenceEntry,
    include_sequence = False
) -> gRPCSequenceEntry:
    """Serialize a `SequenceEntry` as a gRPC object."""
    sequence = ''
    if include_sequence:
        sequence = str(sequence_entry.sequence)

    entry = gRPCSequenceEntry(
        identifier=str(sequence_entry.id),
        is_composite=sequence_entry.is_composite,
        sequence=sequence,
    )
    entry.roles.extend([
        role_dto(role)
        for role in sequence_entry.roles
    ])

    if sequence_entry.parent_links:
        entry.extend([
            entry_link_dto(entry_link)
            for entry_link in sequence_entry.parent_links
        ])

    return entry
