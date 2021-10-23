"""Tests for the JSON output module."""
from unittest import TestCase
import json
from io import StringIO
from Bio.Seq import Seq

from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.plugins.builtin.output.json_output import JSONOutput

class JSONOutputTestCase(TestCase):
    """Test the behavior of SynBioHub provider"""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

        self.compiler_settings = {}
        self.output = JSONOutput(self.fixtures.mk_function_config_object(self.compiler_settings))

    def test_get_assemblies_from_mocked_detail(self):
        """Check that the assemblies serialized by JSON match to original assemblies."""
        assembly_result = self.fixtures.mk_assembly_result_set()

        parts_by_assembly = {}
        sequence_by_assembly = {}
        for assembly in assembly_result:
            parts_by_assembly[assembly.identifier] = [
                part.identifier
                for part in assembly.parts
            ]
            sequence_by_assembly[assembly.identifier] = assembly.sequence

        output_fp = StringIO()
        self.output.output(assembly_result, output_fp)

        result = json.loads(output_fp.getvalue())
        for assembly_result in result['assemblies']:
            self.assertEqual(
                assembly_result['parts'],
                parts_by_assembly[assembly_result['identifier']])

            self.assertEqual(
                assembly_result['sequence'],
                sequence_by_assembly[assembly_result['identifier']])

    def test_check_part_serialization(self):
        """Make sure that the json output of parts matches fixture data."""
        assembly_result = self.fixtures.mk_assembly_result_set()
        expected_parts = set()
        for assembly in assembly_result:
            expected_parts.update([
                part.identifier
                for part in assembly.parts
            ])

        output_fp = StringIO()
        self.output.output(assembly_result, output_fp)
        result = json.loads(output_fp.getvalue())

        result_part_names = set([
            part['identifier']
            for part in result['parts']
        ])
        self.assertEqual(result_part_names, expected_parts)
