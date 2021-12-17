from typing import List
from supergsl.core.sequence import SequenceEntry, Role, EntryLink


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
