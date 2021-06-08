import unittest
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.types.part import Part
from supergsl.core.constants import SO_PROMOTER

class PartTestCase(unittest.TestCase):
    """Test case for SeqPosition."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

    def test_part_roles_are_recapitulated(self):
        _, p1 = self.fixtures.mk_part('SWEET', 300, roles=[
            SO_PROMOTER
        ])

        self.assertEqual(p1.roles, [SO_PROMOTER])
