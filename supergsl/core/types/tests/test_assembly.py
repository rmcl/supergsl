from unittest import TestCase
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures

from supergsl.core.types.assembly import (
    Assembly,
    AssemblyDeclaration,
    AssemblyResultSet
)


class AssemblyResultSetTestCase(TestCase):
    """Test the AssemblyResultSet functionality."""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

    def test_add_assembly(self):
        """Build up a result set by adding assemblies."""
        asm1 = self.fixtures.mk_assembly(identifier='asm1')
        asm2 = self.fixtures.mk_assembly(identifier='asm2')

        result = AssemblyResultSet([])
        result.add_assembly(asm1)
        result.add_assembly(asm2)

        self.assertEqual(list(result), [asm1, asm2])


class AssemblyTestCase(TestCase):
    """Test case for SeqPosition."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()


    def test_assembly_declaration_makes_factors_and_levels(self):
        """Confirm that declared assembly has correct factors and levels."""
        declaration = self.fixtures.mk_assembly_declaration_ex1()

        self.assertEqual(declaration.num_designs, 1*1*1*3)
        factors = [
            (factor.factor_type, [
                str(level)
                for level in factor.levels
            ])
            for factor in declaration.get_factors()
        ]
        self.assertEqual(factors, [
            ('Part', ['uHO']),
            ('Part', ['pGAL0', 'pGAL1', 'pGAL2']),
            ('Part', ['gGENE']),
            ('Part', ['dHO'])
        ])

    def test_assembly_declaration_makes_full_factorial_designs(self):
        """Get the full factorial of designs from a assembly."""
        declaration = self.fixtures.mk_assembly_declaration_ex1()

        part_levels = declaration.get_levels_by_factor_type('Part')
        parts = {
            part.identifier: part
            for part in part_levels
        }

        designs = declaration.get_designs()
        self.assertEqual(list(designs), [
            [parts['uHO'], parts['pGAL0'], parts['gGENE'], parts['dHO']],
            [parts['uHO'], parts['pGAL1'], parts['gGENE'], parts['dHO']],
            [parts['uHO'], parts['pGAL2'], parts['gGENE'], parts['dHO']]
        ])
