from unittest import TestCase
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.constants import (
    SO_PROMOTER,
    PART_SLICE_POSTFIX_START,
    PART_SLICE_POSTFIX_END
)


class PartTestCase(TestCase):
    """Test case for SeqPosition."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

    def test_part_roles_are_recapitulated(self):
        _, p1 = self.fixtures.mk_part('SWEET', 300, roles=[
            SO_PROMOTER
        ])

        self.assertEqual(p1.roles, [SO_PROMOTER])

    def test_part_print(self):
        """Test the part print statement outputs in a nice format."""
        _, part = self.fixtures.mk_part(identifier='BOOM', part_seq_len=25)

        self.assertEqual(
            part.print(),
            "BOOM: %s" % part.sequence)
