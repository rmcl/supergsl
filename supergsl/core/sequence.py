from typing import List, Dict, Optional, Tuple, NamedTuple
from collections import namedtuple
from uuid import UUID, uuid4
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature

from supergsl.core.exception import SequenceStoreError
from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.slice import Slice, Position


class SliceMapping(NamedTuple):
    """Represent the mapping of a sub-slice of a parent sequence onto a new target sequence."""

    parent_entry: 'SequenceEntry'
    source_slice: Slice
    target_slice: Slice


class EntryLink:
    """A link between two sequences in the store or a sequence annotation."""

    def __init__(
        self,
        parent_entry : 'SequenceEntry',
        source_slice : Slice,
        target_slice : Slice,
        annotations
    ):
        self.parent_entry = parent_entry
        self.source = source_slice
        self.target = target_slice
        self.annotations = annotations

    def get_source_slice_sequence(self) -> Seq:
        """Return the sub-sequence of the parent sequence covered by the source slice."""
        sequence_length = self.parent_entry.sequence_length
        start_abs_position = self.source.start.compute_absolute_position(sequence_length)
        end_abs_position = self.source.end.compute_absolute_position(sequence_length)

        if start_abs_position.is_out_of_bounds or end_abs_position.is_out_of_bounds:
            raise SequenceStoreError('OUT OF BOUNDS SEQUENCES NOT CURRENTLY SUPPORTED!')

        parent_sequence = self.parent_entry.sequence
        if start_abs_position.index > end_abs_position.index:
            # This sequence is on the reverse strand. Retrieve the end to start
            # sequence on the forward strand and then take the reverse complement.
            forward_sequence = parent_sequence[end_abs_position.index:start_abs_position.index]
            sequence = forward_sequence.reverse_complement()
        else:
            sequence = parent_sequence[start_abs_position.index:end_abs_position.index]

        return sequence

class SequenceEntry:
    """Represent a sequence in the sequence store."""
    def __init__(
        self,
        store : 'SequenceStore',
        sequence_id: UUID,
        parent_links : Optional[List[EntryLink]] = None,
        reference : Optional[Seq] = None
    ):
        self.sequence_id = sequence_id
        self.parent_links = parent_links
        self.reference = reference

        if not (self.parent_links or self.reference):
            raise SequenceStoreError('Must specify either parents or reference.')
        if self.parent_links and self.reference:
            raise SequenceStoreError('Must only specify parents or a reference sequence, but not both')

    def add_annotation(self, slice, details):
        # Todo: maybe add this function to create EntryLinks for the provided annotations.
        raise NotImplemented()

    @property
    def id(self):
        return self.sequence_id

    @property
    def is_composite(self):
        """SequenceEntries can either by Composite or Hierarchical

        Composite parts are those that derive from more that one part at any point in their sequence lineage.
        Hierarchical parts have only a single parent at every level of their lineage.

        A composite part has been created by the concatination of several reference sequences whereas
        a Hierarhical entry has been sliced out of a *single* reference.

        The primary reason for diferentiating is that it is obvious how to exceed the slice bounds of a hierarchical part.
        For example the part pGAL3[-500:0] would return 500bp upstream of the pGAL3 part's start position.

        In contrast a hierarchical part for example:
            let part = usHO ; pGAL3 ; GENE ; tSDH1 ; dsHO
            part[-500:0]

        It is unclear what the sequence of 500 bp upstream should be? Does the user just want 500 bp of HO or something else?

        """
        if self.reference:
            return False
        if len(self.parent_links) > 1:
            return True

        return self.parent_links[0].is_composite()

    @property
    def reference_sequence(self) -> Seq:
        """Return the `Seq` representing the reference sequence for this position."""
        if self.reference:
            return self.reference

        if self.parent_links:
            if len(self.parent_links) == 1:
                return self.parent_links[0].parent_entry.reference_sequence

        raise SequenceStoreError('Composite parts have multiple reference sequences.')

    @property
    def sequence_length(self) -> int:
        """Return the length of the stored sequence.

        This method is here so we can in the future supply a better implementation
        for figuring this out rather than evaluating the entire sequence entry tree.
        """
        return len(self.sequence)

    @property
    def sequence(self) -> Seq:
        if self.reference:
            return self.reference

        # build_target_sequence_coordinate_map(self.parents)

        sequence_index = 0
        result_sequence = ''

        if not self.parent_links:
            return Seq(result_sequence)

        # TODO: make sure parents are ordered
        for parent_link in self.parent_links:
            start_pos = parent_link.target.start.index
            source_sequence = parent_link.get_source_slice_sequence()

            if start_pos == sequence_index:
                result_sequence += source_sequence
                sequence_index += len(source_sequence)

            elif start_pos < sequence_index:
                overlap = sequence_index - start_pos
                remaining = len(source_sequence) - overlap

                result_sequence += source_sequence[overlap:]
                sequence_index += remaining

            else:
                raise SequenceStoreError(
                    'Gap detected in sequence of SequenceEntry. Entry sequence '
                    'cannot be materialized.')

        return Seq(result_sequence)


    def _build_target_sequence_coordinate_map(self, parent_links):
        """Translate all the target coordinates to be absolulte positions from the five prime start of the Watson strand."""
        for parent_link in parent_links:
            parent_link.target.start
            parent_link.target.end

    def __repr__(self):
        return '%s: %s or %s' % (self.id, self.parent_links, self.reference)


class SequenceStore:

    def __init__(self):
        self._sequences_by_uuid : Dict[UUID, SequenceEntry] = {}

    def lookup(self, sequence_id : UUID) -> SequenceEntry:
        return self._sequences_by_uuid[sequence_id]

    def _create_record_id(self):
        return uuid4()

    def list(self):
        """List sequences in the store."""
        raise NotImplementedError('implement me!')

    def add_from_reference(self, sequence : Seq):

        # TODO: Do we want to support adding from a SeqRecord
        # If yes, consider iterating over the SeqRecord's features and adding them
        # as annotations on EntryLinks

        entry = SequenceEntry(
            store=self,
            sequence_id=self._create_record_id(),
            reference=sequence)

        # Todo: Figure out a better way to do this
        # seq_hash = hash(str(sequence))
        # self._sequences_by_hash[seq_hash]

        self._sequences_by_uuid[entry.id] = entry
        return entry

    def slice(self, sequence_entry : SequenceEntry, sequence_slice : Slice):
        """Create a new entry from a subsequence of the provided sequence_entry.

        The new sequence will be defined by the `sequence_slice`.
        """
        link = EntryLink(
            parent_entry=sequence_entry,
            source_slice=sequence_slice,
            target_slice=Slice(
                Position(0),
                Position(0, relative_to=THREE_PRIME)
            ),
            annotations=None)

        entry = SequenceEntry(
            store=self,
            sequence_id=self._create_record_id(),
            parent_links=[link]
        )
        self._sequences_by_uuid[entry.id] = entry

        return entry

    def concatenate(self, slice_mappings : List[SliceMapping]):
        """Concatenate several sequence slices into a composite sequence."""
        links = []
        for slice_mapping in slice_mappings:
            link = EntryLink(
                slice_mapping.parent_entry,
                slice_mapping.source_slice,
                slice_mapping.target_slice,
                None)
            links.append(link)

        entry_id = self._create_record_id()
        entry = SequenceEntry(self, entry_id, links)
        self._sequences_by_uuid[entry_id] = entry
        return entry
