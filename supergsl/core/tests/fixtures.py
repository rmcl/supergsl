"""Define fixture data for testing SuperGSL."""
from typing import Tuple, List, Optional
from unittest.mock import Mock
from random import randint, choice
from Bio import Restriction
from Bio.Seq import Seq

from supergsl.core.sequence import (
    SequenceStore,
    SequenceEntry,
    SequenceAnnotation,
    SliceMapping,
    Role
)
from supergsl.core.provider import ProviderConfig
from supergsl.core.types.primer import PrimerPair
from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.part import Part
from supergsl.core.types.position import Slice, Position
from supergsl.core.types.builtin import Collection
from supergsl.core.types.assembly import (
    Assembly,
    AssemblyDeclaration,
    AssemblyLevelDeclaration,
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
        symbol_table.insert('uHO', self.mk_part('uHO', 100))

        nested_symbol_table = symbol_table.enter_nested_scope('awesome')
        nested_symbol_table.insert('tHUG', self.mk_part('tHUG', 33))

        return symbol_table

    def mk_random_dna_sequence(
        self,
        seq_len : int,
        excluded_restriction_sites : Optional[List[str]] = None
    ) -> Seq:
        """Make a `Seq` with random DNA of a given length devoid of the given restriction sites."""

        restriction_sites = []
        if excluded_restriction_sites:
            for restriction_site_name in excluded_restriction_sites:
                try:
                    restriction_sites.append(
                        getattr(Restriction, restriction_site_name))
                except AttributeError as error:
                    raise Exception(f'Unknown restriction {restriction_site_name}') \
                        from error

        while True:
            random_seq = Seq(''.join(choice('CGTA') for _ in range(seq_len)))

            bad_sequence = False
            for restriction in restriction_sites:
                if restriction.search(random_seq):
                    bad_sequence = True

            # If we encounter a bad restriction site, go to the top of the loop
            # and generate a new sequence.
            if bad_sequence:
                continue

            # Sequence is free of restriction sites so return.
            return random_seq

    def mk_random_dna_sequence_entry(
        self,
        seq_len : int, annotations : Optional[List[SequenceAnnotation]] = None
    ) -> SequenceEntry:
        """Make a sequence entry with random DNA of a given length."""
        sequence = self.mk_random_dna_sequence(seq_len)
        return self.mk_sequence_entry(sequence, annotations)

    def mk_random_sequence_annotation(self, seq_len : int, roles : List[Role]) -> SequenceAnnotation:
        start = randint(0, int(seq_len / 2))
        end = randint(start, seq_len)
        return SequenceAnnotation.from_five_prime_indexes(
            start,
            end,
            roles=roles,
            payload={
                'pay': 'load'
            })


    def mk_extraction_primers(self, part) -> PrimerPair:
        """Primers are often complementary sequences flanking the DNA sequence of interest.

        Real primers woud use some smarter algorithm to make sure the resulting
        DNA sequence has a melting temperature that is conducive to the parameters
        of a the particular PCR thermocycler.
        """
        print(part.sequence_entry.sequence_length)
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

    def mk_provider_config(self, options=None):
        """Make a config object for instantiating part providers."""
        if not options:
            options = {}
        return ProviderConfig(self.sequence_store, options)

    @property
    def sequence_store(self) -> SequenceStore:
        if not hasattr(self, '_store'):
            self._store = SequenceStore()
        return self._store


    def mk_sequence_entry(
        self,
        sequence : Seq,
        annotations : Optional[List[SequenceAnnotation]] = None
    ) -> SequenceEntry:
        """Create a sequence entry with the given sequence and number of annotaitons."""
        return self.sequence_store.add_from_reference(
            sequence,
            annotations=annotations)

    def mk_sequence_role(self, uri):
        return Role(uri=uri, name=uri, description=uri)

    def mk_part_by_sequence_entry(
        self,
        identifier : str,
        sequence_entry : SequenceEntry,
        mk_primers : bool = True,
        roles : Optional[List[Role]] = None
    ) -> Part:
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

        return part


    def mk_part(
        self,
        identifier : str,
        part_seq_len : int,
        mk_primers : bool = True,
        roles : Optional[List[Role]] = None,
        excluded_restriction_sites : List[str] = None
    ) -> Part:
        """Create a mock Part.

        Part is derived from a reference sequence three times longer than the part
        itself. That parent sequence is sliced and the part returned is the slice:
        100:100+part_seq_len.
        """

        # Ensure that the reference sequence length is at least
        # 100bp + part_seq_len. That way when we slice the part out 100bp in we
        # will avoid an out of bounds error.
        ref_seq_len = max(part_seq_len * 3, 100 + part_seq_len)
        reference_sequence = self.mk_random_dna_sequence(
            ref_seq_len,
            excluded_restriction_sites=excluded_restriction_sites)

        sequence_entry = self.mk_sequence_entry(reference_sequence)
        part_sequence_entry = self.sequence_store.slice(
            sequence_entry,
            Slice.from_five_prime_indexes(100, 100 + part_seq_len))

        return self.mk_part_by_sequence_entry(
            identifier,
            part_sequence_entry,
            mk_primers,
            roles)


    def mk_part_collection(self, num_parts=3, part_len=100) -> Collection:
        """Return a `Collection` of parts."""

        return Collection([
            self.mk_part('pGAL%d' % index, part_len)
            for index in range(num_parts)
        ])

    def mk_assembly_declaration_ex1(self, name='test_declaration') -> AssemblyDeclaration:
        """Make a AssemblyDeclaration with 4 parts including one with a collection of three parts"""
        promoter_collection = self.mk_part_collection(num_parts=3)
        gene = self.mk_part('gGENE', 500)
        upstream = self.mk_part('uHO', 50)
        downstream = self.mk_part('dHO', 50)

        return AssemblyDeclaration(name, [
            AssemblyLevelDeclaration(upstream, None),
            AssemblyLevelDeclaration(promoter_collection, None),
            AssemblyLevelDeclaration(gene, None),
            AssemblyLevelDeclaration(downstream, None)
        ])

    def mk_assembly(self, identifier='asm1', num_parts=2) -> Assembly:
        """Create a `Assembly` containing num_parts with random sequences of len 100 to 1000."""
        parts : List[Part] = list([
            self.mk_part('part-%03d' % part_index, randint(100, 1000))
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
                num_parts=randint(1,5))
            for assembly_index in range(num_assembly)
        ])
