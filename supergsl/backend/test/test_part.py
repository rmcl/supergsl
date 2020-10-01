import mock
from unittest import TestCase
from supergsl.core.test.fixtures import SuperGSLCoreFixtures
from supergsl.backend.parts import (
    SliceAndBuildPartSequencePass,
    DNASlice,
    SequencePosition,
    Part
)


class PartTestCase(TestCase):

    @property
    def fixtures(self):
        return SuperGSLCoreFixtures()

    def test_build_part_type_slice(self):
        p = SliceAndBuildPartSequencePass()
        source_part = mock.Mock()
        dna_slice = p.build_part_type_slice(source_part, 'promoter')

        self.assertEquals(dna_slice.source_part, source_part)
        self.assertEquals(dna_slice.part_type, 'promoter')
        self.assertEquals(dna_slice.left, SequencePosition(
            x=500, rel_to='FivePrime', approximate=True)
        )
        self.assertEquals(dna_slice.right, SequencePosition(
            x=-1, rel_to='FivePrime', approximate=False)
        )

    def test_slice_part_by_part_type(self):

        p = SliceAndBuildPartSequencePass()
        source_part = Part(
            name='pTEST',
            sequence='',
            parent_part=None
        )

        dna_slice = p.build_part_type_slice(source_part, 'promoter')

        print(dna_slice)
        self.assertEquals(True, False)
