from typing import List
from supergsl.core.sequence import SequenceEntry, SequenceAnnotation, Role, EntryLink


def filter_links_by_roles(sequence_entry : SequenceEntry, roles : List[Role]) -> List[EntryLink]:
    """Filter Sequence Entry links by their roles."""

    if not sequence_entry.parent_links:
        return []

    included_links : List[EntryLink] = []
    for parent_link in sequence_entry.parent_links:
        include = False
        for role in parent_link.roles:
            if role in roles:
                include = True

        if include:
            included_links.append(parent_link)

    return included_links

def filter_annotations_by_roles(
    sequence_entry : SequenceEntry,
    roles : List[Role]
) -> List[SequenceAnnotation]:
    """Return annotaitons for a sequence entry with that have one of the specified roles."""
    included_annotations = []
    for annotation in sequence_entry.annotations():
        include = False
        if not annotation.roles:
            continue

        for annotation_role in annotation.roles:
            if annotation_role in roles:
                include = True

        if include:
            included_annotations.append(annotation)

    return included_annotations
