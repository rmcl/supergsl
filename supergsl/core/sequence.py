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
from supergsl.core.types.position import Slice, Position, AbsoluteSlice, AbsolutePosition


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

    def __eq__(self, other):
        """Roles with the same uri should be considered identical."""
        return self.uri == other.uri


class SliceMapping(NamedTuple):
    """Represent the mapping of a sub-slice of a parent sequence onto a new target sequence."""
    parent_entry: 'SequenceEntry'
    source_slice: Slice
    target_slice: Slice
    roles : Optional[List[Role]] = None


class SequenceAnnotation(NamedTuple):
    """Store a annotation at a position relative to a sequence."""
    location: Slice
    roles: Optional[List[Role]]
    payload: dict

    def __eq__(self, other):
        """Annotations with the same position, roles and payload are the same."""
        return (
            self.location == other.location and
            self.roles == other.roles and
            self.payload == other.payload
        )

    def derive_absolute_position_annotation(
        self,
        annotation_start : AbsolutePosition
    ) -> 'SequenceAnnotation':
        """Derive a new Annotation with a location relative to the given start position."""

        # TODO: Highly skeptical that this will work on the reverse strand.
        # Need to build test case.
        annotation_location = Slice(
            Position(
                self.location.start.index - annotation_start.index,
                self.location.start.relative_to,
                self.location.start.approximate),
            Position(
                self.location.end.index - annotation_start.index,
                self.location.end.relative_to,
                self.location.end.approximate)
        )

        print('ARGH', annotation_location, annotation_start.index)
        return SequenceAnnotation(
            location=annotation_location,
            roles=self.roles,
            payload=self.payload)

    @classmethod
    def from_five_prime_indexes(
        cls,
        start : int,
        end : int,
        roles : Optional[List[Role]],
        payload : dict
    ):
        """Create a sequence annotation relative to five prime end of a sequence."""
        return SequenceAnnotation(
            Slice.from_five_prime_indexes(start, end),
            roles,
            payload
        )


class EntryLink:
    """A link between two sequences in the store."""

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

    def sequence_annotations(self) -> List[SequenceAnnotation]:
        """Return all annotations within the slice of the parent entry."""
        return self.parent_entry.sequence_annotations_for_slice(self.source_slice)


class SequenceEntry:
    """Represent a sequence in the sequence store."""
    def __init__(
        self,
        sequence_store : 'SequenceStore',
        sequence_id: UUID,
        roles : Optional[List[Role]] = None,
        parent_links : Optional[List[EntryLink]] = None,
        reference : Optional[Seq] = None,
        external_ids : Optional[Dict[str, str]] = None,
        annotations : Optional[List[SequenceAnnotation]] = None
    ):
        self.sequence_store = sequence_store
        self.sequence_id = sequence_id
        self.roles = roles if roles else []
        self.parent_links = parent_links
        self.reference = reference
        self._annotations = annotations or []

        if not external_ids:
            self._external_ids : Dict[str, str] = {}

        if not (self.parent_links or self.reference):
            raise SequenceStoreError('Must specify either parents or reference.')
        if self.parent_links and self.reference:
            raise SequenceStoreError(
                'Must only specify parents or a reference sequence, but not both')

    def _get_local_sequence_annotations_for_slice(self, desired_slice : Slice):
        """Get annotations for a slice of this entity."""

        # TODO MAYBE WE DONT WANT TO CONVERT TO ABSOLUTE SLICE HERE!
        # this is going to be called through a entity link so it would be better
        # if this returned relative positions that were then converted.
        absolute_desired_slice = desired_slice.build_absolute_slice(self.sequence_length)

        filtered_annotations = []
        for annotation in self._annotations:
            annotation_absolute_slice = annotation.location.build_absolute_slice(
                self.sequence_length)

            ### TODO: VERY SKEPTICAL OF THIS LOGIC!
            ### REVISIT AND MAKE SURE IT DOES WHAT WE EXPECT.
            if annotation_absolute_slice.start < absolute_desired_slice.start:
                continue
            if annotation_absolute_slice.end >= absolute_desired_slice.end:
                continue

            filtered_annotations.append(annotation)

        return filtered_annotations

    def sequence_annotations_for_slice(self, desired_slice: Slice) -> List[SequenceAnnotation]:
        """Retrieve annotations contained in the given slice."""

        annotations = []
        # Retrieve annotations stored on this entry
        local_annotations = self._get_local_sequence_annotations_for_slice(desired_slice)
        for annotation in local_annotations:
            absolute_slice = annotation.location.build_absolute_slice(self.sequence_length)
            annotations.append(annotation)

        if self.parent_links:
            for parent_link in self.parent_links:
                """
                absolute_parent_target_slice = parent_link.target_slice.build_absolute_slice(
                    self.sequence_length)
                """

                target_start_pos = AbsolutePosition(
                    self.sequence_length,
                    parent_link.source_slice.start.index,
                    False)

                for parent_annotation in parent_link.sequence_annotations():
                    print('YOOO111', self.sequence_length)

                    new_annotation = parent_annotation.derive_absolute_position_annotation(target_start_pos)

                    ## PROBLEM IS WE NEED TO DERIVE THIS LOCATION RELATIVE TO THE START
                    # OF THE PARENT PART IN THE CHILD PART
                    print(parent_annotation.location, new_annotation.location)

                annotations.append(new_annotation)

        annotations = sorted(
            annotations,
            key=lambda annotation: annotation.location.start.index)

        return annotations

    def sequence_annotations(self) -> List[SequenceAnnotation]:
        """Return the sequence annotations for this entry.

        Recursively include all annotations present in the parent links composing this entry.
        """
        return self.sequence_annotations_for_slice(Slice.from_entire_sequence())

    def add_annotation(self, annotation : SequenceAnnotation):
        """Add an annotation to the entry."""
        self._annotations.append(annotation)

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
        sequence : Union[Seq, SeqRecord],
        roles : Optional[List[Role]] = None,
        external_ids : Optional[Dict[str,str]] = None,
        annotations : Optional[List[SequenceAnnotation]] = None
    ):
        """Add a sequence to the store."""

        # TODO: Do we want to support adding from a SeqRecord
        # If yes, consider iterating over the SeqRecord's features and adding them
        # as annotations on EntryLinks
        if isinstance(sequence, SeqRecord):
            sequence_record = sequence
            sequence = sequence_record.seq

        external_ids = external_ids or {}
        self._check_for_duplicate_external_ids(external_ids)

        entry = SequenceEntry(
            sequence_store=self,
            sequence_id=self.__create_record_id(),
            reference=sequence,
            roles=roles if roles else [],
            external_ids=external_ids,
            annotations=annotations)

        # Index any provided external ids.
        self._add_entry_to_external_id_index(entry, external_ids)

        # Index the entry's UUID
        self._sequences_by_uuid[entry.id] = entry
        return entry

    def slice(
        self,
        sequence_entry : SequenceEntry,
        sequence_slice : Slice,
        new_sequence_roles : Optional[List[Role]] = None,
        entry_link_roles : Optional[List[Role]] = None,
        annotations : Optional[List[SequenceAnnotation]] = None
    ) -> SequenceEntry:
        """Create a new entry from a subsequence of the provided sequence_entry.

        The new sequence will be defined by the `sequence_slice`.
        """
        new_sequence_roles = new_sequence_roles or []
        entry_link_roles = entry_link_roles or []

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
            roles=new_sequence_roles,
            annotations=annotations
        )
        self._sequences_by_uuid[entry.id] = entry

        return entry

    def slice_from_annotation(
        self,
        sequence_entry : SequenceEntry,
        annotation : SequenceAnnotation
    ) -> SequenceEntry:
        """Create a new sequence entry for an annotation."""
        return self.slice(
            sequence_entry,
            annotation.location,
            new_sequence_roles=annotation.roles)

    def _check_for_duplicate_external_ids(self, external_ids):
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

    def _add_entry_to_external_id_index(self, entry : SequenceEntry, external_ids : Dict[str,str]):
        """Index any provided external ids."""
        for external_system_name, external_id in external_ids.items():
            if external_system_name not in self._sequences_by_external_id:
                self._sequences_by_external_id[external_system_name] = {}

            self._sequences_by_external_id[external_system_name][external_id] = entry

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

        if not new_sequence_roles:
            new_sequence_roles = []

        entry_id = self.__create_record_id()
        entry = SequenceEntry(self, entry_id, new_sequence_roles, links)
        self._sequences_by_uuid[entry_id] = entry
        return entry
