from typing import List, Dict, Optional, Tuple, NamedTuple, Union
from uuid import UUID, uuid4
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from supergsl.core.exception import (
    SequenceStoreError,
    SequenceNotFoundError,
    DuplicateSequenceError
)
from supergsl.core.constants import THREE_PRIME, STRAND_CRICK
from supergsl.core.types.position import Slice, Position, AbsoluteSlice


def get_slice_sequence_from_reference(sequence_reference, absolute_slice) -> Seq:
    """Return the sub-sequence of the reference sequence covered by the provided slice."""

    if absolute_slice.start.is_out_of_bounds or absolute_slice.end.is_out_of_bounds:
        raise SequenceStoreError(
            'The position defined in the slice exceeds the bounds of the reference sequence.')

    if absolute_slice.strand == STRAND_CRICK:
        # This sequence is on the reverse strand. Retrieve the end to start
        # sequence on the forward strand and then take the reverse complement.
        reverse_sequence = sequence_reference.reverse_complement()
        sequence = reverse_sequence[absolute_slice.start.index:absolute_slice.end.index]
    else:
        sequence = sequence_reference[absolute_slice.start.index:absolute_slice.end.index]

    return sequence


class Role(NamedTuple):
    """Represent a sequence role."""
    uri : str
    name : str
    description : str


class SliceMapping(NamedTuple):
    """Represent the mapping of a sub-slice of a parent sequence onto a new target sequence."""

    parent_entry: 'SequenceEntry'
    source_slice: Slice
    target_slice: Slice
    roles : Optional[List[Role]] = None


class EntryLink:
    """A link between two sequences in the store or a sequence annotation."""

    def __init__(
        self,
        parent_entry : 'SequenceEntry',
        source_slice : Slice,
        target_slice : Slice,
        roles : List[Role]
    ):
        self.parent_entry = parent_entry
        self.source_slice = source_slice
        self.target_slice = target_slice
        self.roles = roles

    @property
    def sequence(self) -> Seq:
        """Return a sequence representing the slice of the parent entry."""
        reference_source_sequence, absolute_source_slice = \
            self.parent_entry.get_slice_absolute_position_and_reference(
                self.source_slice)

        return get_slice_sequence_from_reference(
            reference_source_sequence,
            absolute_source_slice)


class SequenceEntry:
    """Represent a sequence in the sequence store."""
    def __init__(
        self,
        sequence_store : 'SequenceStore',
        sequence_id: UUID,
        roles : Optional[List[Role]] = None,
        parent_links : Optional[List[EntryLink]] = None,
        reference : Optional[Seq] = None,
        external_ids : Optional[Dict[str, str]] = None
    ):
        self.sequence_store = sequence_store
        self.sequence_id = sequence_id
        self.roles = roles if roles else []
        self.parent_links = parent_links
        self.reference = reference

        if not external_ids:
            self._external_ids : Dict[str, str] = {}

        if not (self.parent_links or self.reference):
            raise SequenceStoreError('Must specify either parents or reference.')
        if self.parent_links and self.reference:
            raise SequenceStoreError(
                'Must only specify parents or a reference sequence, but not both')

    @property
    def external_ids(self):
        """Return a dictionary containing all external ids."""
        return self._external_ids

    def get_external_id(self, external_system_name : str):
        """Return the external id for a specific system."""
        return self._external_ids[external_system_name]

    @property
    def id(self) -> UUID:
        """Return the UUID of the sequence entry."""
        return self.sequence_id

    @property
    def is_composite(self):
        """SequenceEntries can either by Composite or Hierarchical

        Composite parts are those that derive from more that one part at any
        point in their sequence lineage. Hierarchical parts have only a single
        parent at every level of their lineage.

        A composite part has been created by the concatination of several
        reference sequences whereas a Hierarhical entry has been sliced out of
        a *single* reference.

        The primary reason for diferentiating is that it is obvious how to exceed
        the slice bounds of a hierarchical part. For example the part pGAL3[-500:0]
        would return 500bp upstream of the pGAL3 part's start position.

        In contrast a hierarchical part for example:
            let part = usHO ; pGAL3 ; GENE ; tSDH1 ; dsHO
            part[-500:0]

        It is unclear what the sequence of 500 bp upstream should be? Does the
        user just want 500 bp of HO or something else?

        """
        if self.reference:
            return False
        if len(self.parent_links) > 1:
            return True

        return False


    @property
    def sequence_length(self) -> int:
        """Return the length of the stored sequence."""
        entire_sequence_slice = Slice.from_entire_sequence()
        _, absolute_slice = self.get_slice_absolute_position_and_reference(
            entire_sequence_slice)
        return len(absolute_slice)


    def get_slice_absolute_position_and_reference(
        self,
        target_slice : Slice
    ) -> Tuple[Seq, AbsoluteSlice]:
        """Return the reference sequence and absolute position position of slice.

        SOURCE SEQ: 5'----------------3'
        TARGET SEQ: 5'--|XXX|---------3'

        source_slice: <start 5':0, end: 3'0>
        target slice: <start: 5':3, end: 5'6>
        reference: "SOURCE SEQ"
        """

        if self.reference:
            absolute_slice = target_slice.build_absolute_slice(len(self.reference))
            return self.reference, absolute_slice

        if not self.is_composite:
            if not self.parent_links:
                raise Exception('no parent links')

            parent_link = self.parent_links[0]
            parent_sequence_entry = parent_link.parent_entry
            reference_sequence, parent_absolute_slice = \
                parent_sequence_entry.get_slice_absolute_position_and_reference(
                    parent_link.source_slice)

            absolute_source_slice = parent_absolute_slice.derive_from_relative_slice(target_slice)
            return reference_sequence, absolute_source_slice

        # This thing is a composite part so lets just materialize the sequence.
        reference = self.sequence
        absolute_slice = target_slice.build_absolute_slice(len(reference))
        return reference, absolute_slice


    @property
    def sequence(self) -> Seq:
        """Return the complete sequence of this sequence entry."""
        if self.reference:
            return self.reference

        result_sequence = ''
        sequence_index = 0
        if not self.parent_links:
            return Seq(result_sequence)

        # Iterate over parent links and concatenate sequences.
        for parent_link in self.parent_links:

            reference_source_sequence, absolute_source_slice = \
                parent_link.parent_entry.get_slice_absolute_position_and_reference(
                    parent_link.source_slice)

            source_sequence = get_slice_sequence_from_reference(
                reference_source_sequence,
                absolute_source_slice)

            target_start_pos = parent_link.target_slice.start.index

            if target_start_pos == sequence_index:
                result_sequence += str(source_sequence)
                sequence_index += len(source_sequence)

            elif target_start_pos < sequence_index:
                overlap = sequence_index - target_start_pos
                remaining = len(source_sequence) - overlap

                result_sequence += str(source_sequence[overlap:])
                sequence_index += remaining

            else:
                raise SequenceStoreError(
                    'Gap detected in sequence of SequenceEntry. Entry sequence '
                    'cannot be materialized.')

        return Seq(result_sequence)

    def __repr__(self):
        return '%s: %s or %s' % (self.id, self.parent_links, self.reference)


class SequenceStore:
    """Store sequences and provided methods for deriving new sequences from existing ones."""

    def __init__(self):
        self._sequences_by_uuid : Dict[UUID, SequenceEntry] = {}
        self._links_by_parent_entry_id : Dict[UUID, List[EntryLink]] = {}

        """
        Maintain a index of external unique ids.

        Example:
         {
           'uniprot': {
               'Q00955': SequenceEntry,
               ....
            },
            ...
        }
        """
        self._sequences_by_external_id : Dict[str, Dict[str, SequenceEntry]] = {}
        self._links_by_uuid : Dict[UUID, List[EntryLink]] = {}


    def lookup(self, sequence_id : UUID) -> SequenceEntry:
        """Lookup `SequenceEntry` by it's id."""
        return self._sequences_by_uuid[sequence_id]

    def lookup_by_external_id(self, external_system_name : str, external_id : str) -> SequenceEntry:
        """Find sequences with a given external id."""
        try:
            external_system = self._sequences_by_external_id[external_system_name]
        except KeyError as external_sys_error:
            raise SequenceNotFoundError('Unknown external system "%s"' % (
                external_system_name)) from external_sys_error

        try:
            return external_system[external_id]
        except KeyError as external_id_error:
            raise SequenceNotFoundError('Unknown sequence "%s" in "%s"' % (
                external_id, external_system_name)) from external_id_error

    def __create_record_id(self):
        """Create a new random id for a sequence entry."""
        return uuid4()

    def list(self):
        """List sequences in the store."""
        return self._sequences_by_uuid.values()

    def add_from_reference(
        self,
        sequence : Seq,
        roles : Optional[List[Role]] = None,
        external_ids : Optional[Dict[str,str]] = None):
        """Add a sequence to the store."""

        # TODO: Do we want to support adding from a SeqRecord
        # If yes, consider iterating over the SeqRecord's features and adding them
        # as annotations on EntryLinks
        if isinstance(sequence, SeqRecord):
            sequence_record = sequence
            sequence = sequence_record.seq

        if not external_ids:
            external_ids = {}

        for external_system_name, external_id in external_ids.items():
            try:
                self.lookup_by_external_id(external_system_name, external_id)
            except SequenceNotFoundError:
                continue
            else:
                raise DuplicateSequenceError(
                    'Sequence "%s" from "%s" already exists in the store.' % (
                        external_id,
                        external_system_name
                    ))

        entry = SequenceEntry(
            sequence_store=self,
            sequence_id=self.__create_record_id(),
            reference=sequence,
            roles=roles if roles else [],
            external_ids=external_ids)

        # Index any provided external ids.
        for external_system_name, external_id in external_ids.items():
            if external_system_name not in self._sequences_by_external_id:
                self._sequences_by_external_id[external_system_name] = {}

            self._sequences_by_external_id[external_system_name][external_id] = entry

        # Index the entry's UUID
        self._sequences_by_uuid[entry.id] = entry
        return entry

    def slice(
        self,
        sequence_entry : SequenceEntry,
        sequence_slice : Slice,
        new_sequence_roles : Optional[List[Role]] = None,
        entry_link_roles : Optional[List[Role]] = None
    ):
        """Create a new entry from a subsequence of the provided sequence_entry.

        The new sequence will be defined by the `sequence_slice`.
        """
        if not new_sequence_roles:
            new_sequence_roles = []

        if not entry_link_roles:
            entry_link_roles = []

        link = EntryLink(
            parent_entry=sequence_entry,
            source_slice=sequence_slice,
            target_slice=Slice(
                Position(0),
                Position(0, relative_to=THREE_PRIME)
            ),
            roles=entry_link_roles)

        entry = SequenceEntry(
            sequence_store=self,
            sequence_id=self.__create_record_id(),
            parent_links=[link],
            roles=new_sequence_roles
        )
        self._sequences_by_uuid[entry.id] = entry
        self._add_link_to_parent_index(link)

        return entry

    def _add_link_to_parent_index(self, entry_link : EntryLink):
        parent_entry_id = entry_link.parent_entry.id
        if parent_entry_id not in self._links_by_parent_entry_id:
            self._links_by_parent_entry_id[parent_entry_id] = []
        self._links_by_parent_entry_id[parent_entry_id].append(entry_link)

    def concatenate(
        self,
        slice_mappings : List[SliceMapping],
        new_sequence_roles : Optional[List[Role]] = None
    ) -> SequenceEntry:
        """Concatenate several sequence slices into a composite sequence."""
        links = []
        for slice_mapping in slice_mappings:
            link = EntryLink(
                slice_mapping.parent_entry,
                slice_mapping.source_slice,
                slice_mapping.target_slice,
                slice_mapping.roles if slice_mapping.roles else [])
            links.append(link)

            self._add_link_to_parent_index(link)

        if not new_sequence_roles:
            new_sequence_roles = []

        entry_id = self.__create_record_id()
        entry = SequenceEntry(self, entry_id, new_sequence_roles, links)
        self._sequences_by_uuid[entry_id] = entry
        return entry
