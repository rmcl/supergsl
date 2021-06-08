"""Unit tests for the Fuse Assembler."""
import unittest
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.plugins.builtin.fuse import FusionAssembler

class FuseAssemblerTestCase(unittest.TestCase):
    """Test that the parser correctly tokens into a valid AST."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.assembler = FusionAssembler({})

    def test_assemble_concatenates_part_sequences(self):
        """Fuse assembler is expected to append part sequences together."""
        declaration = self.fixtures.mk_assembly_declaration_ex1()
        parts = {
            part.identifier: part
            for part in declaration.get_levels_by_factor_type('Part')
        }

        result = list(self.assembler.assemble([declaration]))

        expected_results = [
            (
                'ASM-test_declaration-000',
                [
                    parts['uHO'],
                    parts['pGAL0'],
                    parts['gGENE'],
                    parts['dHO']
                ]
            ),
            (
                'ASM-test_declaration-001',
                [
                    parts['uHO'],
                    parts['pGAL1'],
                    parts['gGENE'],
                    parts['dHO']
                ]
            ),
            (
                'ASM-test_declaration-002',
                [
                    parts['uHO'],
                    parts['pGAL2'],
                    parts['gGENE'],
                    parts['dHO']
                ]
            )
        ]

        self.assertEqual(len(result), 3)
        for result_index, expected_result in enumerate(expected_results):
            expected_identifier = expected_result[0]
            expected_parts = expected_result[1]
            self.assertEqual(result[result_index].identifier, expected_identifier)

            self.assertEqual(result[result_index].get_required_parts(), expected_parts)
            self.assertEqual(result[result_index].get_sequence(), ''.join([
                str(part.get_sequence().seq)
                for part in expected_parts
            ]))
