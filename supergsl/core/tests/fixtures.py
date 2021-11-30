"""Define fixture data for testing SuperGSL."""
from typing import Tuple, List, Optional
from unittest.mock import Mock
import random
from Bio.Seq import Seq
from supergsl.core.sequence import (
    SequenceStore,
    SequenceEntry,
    SliceMapping,
    Role
)
from supergsl.core.function import SuperGSLFunctionConfig
from supergsl.core.parts.provider import PartProviderConfig
from supergsl.core.types.primer import PrimerPair
from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.part import Part
from supergsl.core.types.position import Slice, Position
from supergsl.core.types.builtin import Collection
from supergsl.core.types.assembly import (
    Assembly,
    AssemblyDeclaration,
    AssemblyResultSet
)
from supergsl.core.symbol_table import SymbolTable

class SuperGSLCoreFixtures(object):

    def mk_global_symbol_table(self) -> SymbolTable:
        """Instantiate a symbol table with a sequencestore."""
        symbol_table = SymbolTable('global', None)
        symbol_table.insert('sequences', self.sequence_store)
        return symbol_table

    def mk_symbol_table(self) -> SymbolTable:
        """Create a simple symbol table with a nested scope."""
        symbol_table = SymbolTable('awesome', None)
        symbol_table.insert('uHO', self.mk_part('uHO', 100)[1])

        nested_symbol_table = symbol_table.enter_nested_scope('awesome')
        nested_symbol_table.insert('tHUG', self.mk_part('tHUG', 33)[1])

        return symbol_table

    def mk_random_dna_sequence(self, seq_len : int) -> Seq:
        """Make a `Seq` with random DNA of a given length."""
        seq_str = ''.join(random.choice('CGTA') for _ in range(seq_len))
        return Seq(seq_str)

    def mk_random_dna_sequence_entry(self, seq_len : int) -> SequenceEntry:
        """Make a sequence entry with random DNA of a given length."""
        sequence = self.mk_random_dna_sequence(seq_len)
        return self.mk_sequence_entry(sequence)

    def mk_extraction_primers(self, part) -> PrimerPair:
        """Primers are often complementary sequences flanking the DNA sequence of interest.

        Real primers woud use some smarter algorithm to make sure the resulting
        DNA sequence has a melting temperature that is conducive to the parameters
        of a the particular PCR thermocycler.
        """
        complement_sequence = part.sequence.complement()
        complement_sequence_entry = self.mk_sequence_entry(complement_sequence)

        # complement_sequence[:20]
        forward = self.sequence_store.slice(
            complement_sequence_entry,
            Slice.from_five_prime_indexes(0,20)
        )
        # complement_sequence[-20:]
        reverse = self.sequence_store.slice(
            complement_sequence_entry,
            Slice(Position(20), Position(0, relative_to=THREE_PRIME)))

        return PrimerPair.from_sequence_entries(forward, reverse)

    def mk_function_config_object(self, options=None):
        """Make a config object for instantiating supergsl functions."""
        if not options:
            options = {}

        return SuperGSLFunctionConfig(self.sequence_store, options)

    def mk_part_provider_config(self, options=None):
        """Make a config object for instantiating part providers."""
        if not options:
            options = {}
        return PartProviderConfig(self.sequence_store, options)

    @property
    def sequence_store(self) -> SequenceStore:
        if not hasattr(self, '_store'):
            self._store = SequenceStore()
        return self._store

    def mk_sequence_entry(self, sequence : Seq) -> SequenceEntry:
        return self.sequence_store.add_from_reference(sequence)

    def mk_sequence_role(self, uri):
        return Role(uri=uri, name=uri, description=uri)

    def mk_part_by_sequence_entry(
        self,
        identifier : str,
        sequence_entry : SequenceEntry,
        mk_primers : bool = True,
        roles : Optional[List[Role]] = None
    ):
        """Create a mock Part from a given SequenceEntry."""
        part = Part(
            identifier,
            sequence_entry=sequence_entry,
            provider=Mock(),
            roles=roles
        )

        if mk_primers:
            primer_pair = self.mk_extraction_primers(part)
            part.set_extraction_primers(primer_pair)

        return sequence_entry.sequence, part


    def mk_part(
        self,
        identifier : str,
        part_seq_len : int,
        mk_primers : bool = True,
        roles : Optional[List[Role]] = None
    ) -> Tuple[Seq,Part]:
        """Create a mock Part.

        Part is derived from a reference sequence three times longer than the part
        itself.
        """
        ref_seq_len = part_seq_len * 3
        reference_sequence = self.mk_random_dna_sequence(ref_seq_len)
        sequence_entry = self.mk_sequence_entry(reference_sequence)
        return self.mk_part_by_sequence_entry(
            identifier,
            sequence_entry,
            mk_primers,
            roles)


    def mk_part_collection(self, num_parts=3, part_len=100) -> Collection:
        """Return a `Collection` of parts."""

        return Collection([
            self.mk_part('pGAL%d' % index, part_len)[1]
            for index in range(num_parts)
        ])

    def mk_assembly_declaration_ex1(self, name='test_declaration') -> AssemblyDeclaration:
        """Make a AssemblyDeclaration with 4 parts including one with a collection of three parts"""
        promoter_collection = self.mk_part_collection(num_parts=3)
        gene = self.mk_part('gGENE', 500)[1]
        upstream = self.mk_part('uHO', 50)[1]
        downstream = self.mk_part('dHO', 50)[1]

        return AssemblyDeclaration(name, [
            upstream,
            promoter_collection,
            gene,
            downstream
        ])

    def mk_assembly(self, identifier='asm1', num_parts=2) -> Assembly:
        """Create a `Assembly` containing num_parts with random sequences of len 0 to 100."""
        parts : List[Part] = list([
            self.mk_part('part-%03d' % part_index, random.randint(20, 100))[1]
            for part_index in range(num_parts)
        ])

        cur_seq_pos = 0
        slice_mappings : List[SliceMapping] = []
        part_mappings : List[Tuple[Part, Slice]] = []

        for part in parts:
            start_pos = cur_seq_pos
            end_pos = cur_seq_pos + len(part.sequence)
            cur_seq_pos = end_pos

            source_slice = Slice.from_entire_sequence()
            target_slice = Slice.from_five_prime_indexes(start_pos, end_pos)

            slice_mappings.append(
                SliceMapping(part.sequence_entry, source_slice, target_slice))
            part_mappings.append(
                (part, target_slice))

        assembly_sequence_entry = self.sequence_store.concatenate(slice_mappings)

        new_part = Part(
            identifier,
            assembly_sequence_entry,
            provider=self,
            description='This is a great assembly named %s' % identifier)

        assembly = Assembly(
            identifier,
            new_part,
            reagents=parts,
            description='This is a great assembly named %s' % identifier)

        return assembly

    def mk_assembly_result_set(self, num_assembly=2) -> AssemblyResultSet:
        """Create an AssemblyResultSet with `num_assembly` assemblies."""
        return AssemblyResultSet([
            self.mk_assembly(
                'asm-%d' % assembly_index,
                num_parts=random.randint(1,5))
            for assembly_index in range(num_assembly)
        ])
