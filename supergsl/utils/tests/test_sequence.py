import sys
import textwrap
from io import StringIO
from unittest import TestCase

from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.sequence import SliceMapping
from supergsl.core.types.position import Slice
from supergsl.utils.sequence import (
    filter_links_by_roles,
    filter_annotations_by_roles
)


class SequenceUtilitiesTestCases(TestCase):
    """Test sequence related helper methods."""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.sequence_store = self.fixtures.sequence_store

    def test_filter_links_by_roles(self):
        """Test that we can filter parent links by a list of roles."""
        role1 = self.fixtures.mk_sequence_role('AWESOME')
        role2 = self.fixtures.mk_sequence_role('AN_OKAY_ROLE')

        entry1 = self.fixtures.mk_random_dna_sequence_entry(200)
        entry2 = self.fixtures.mk_random_dna_sequence_entry(150)

        new_entry = self.sequence_store.concatenate([
            SliceMapping(
                parent_entry=entry1,
                source_slice=Slice.from_entire_sequence(),
                target_slice=Slice.from_five_prime_indexes(0, 200),
                roles=[role1]),
            SliceMapping(
                parent_entry=entry2,
                source_slice=Slice.from_entire_sequence(),
                target_slice=Slice.from_five_prime_indexes(200, 200+150),
                roles=[role1, role2])
        ])

        results = filter_links_by_roles(new_entry, [])
        self.assertEqual(len(results), 0)

        results_2 = filter_links_by_roles(new_entry, [role1])
        self.assertEqual(len(results_2), 2)
        self.assertEqual({
            result.parent_entry
            for result in results_2
        }, { entry1, entry2 })

        results_3 = filter_links_by_roles(new_entry, [role2])
        self.assertEqual(len(results_3), 1)
        self.assertEqual({
            result.parent_entry
            for result in results_3
        }, { entry2 })

    def test_filter_annotations_by_roles(self):
        """Test that we can filter annotations on a sequence entry by roles."""
        role1 = self.fixtures.mk_sequence_role('AWESOME')
        role2 = self.fixtures.mk_sequence_role('AN_OKAY_ROLE')
        role3 = self.fixtures.mk_sequence_role('BOOP_ROLE')

        annotations = [
            self.fixtures.mk_random_sequence_annotation(200, roles=[role1]),
            self.fixtures.mk_random_sequence_annotation(200, roles=[role1, role2]),
            self.fixtures.mk_random_sequence_annotation(200, roles=[role3])
        ]
        entry1 = self.fixtures.mk_random_dna_sequence_entry(200, annotations)

        self.assertEqual(
            filter_annotations_by_roles(entry1, [role2]),
            [annotations[1]])

        self.assertEqual(
            filter_annotations_by_roles(entry1, [role3]),
            [annotations[2]])

        self.assertEqual(
            len(filter_annotations_by_roles(entry1, [role1, role3])),
            3)
