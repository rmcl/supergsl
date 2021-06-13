"""Tests for the SBOL output module."""
import unittest
import requests
from io import StringIO
from unittest.mock import Mock, call, patch
from Bio.Seq import Seq

from sbol2 import Document

from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.plugins.builtin.output.sbol_output import SBOLOutput


class SBOLOutputTestCase(unittest.TestCase):
    """Test the behavior of SynBioHub provider"""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.sbol_output = SBOLOutput({})

    @patch('supergsl.plugins.builtin.output.sbol_output.Document')
    def test_output(self, sbol_document_class_mock):
        """Test that output method makes appropriate calls to setup SBOL doc and add assemblies."""
        sbol_doc_mock = sbol_document_class_mock.return_value
        sbol_doc_mock.writeString.return_value = 'THE SBOL DOC!'

        results = self.fixtures.mk_assembly_result_set()
        assembly_fixtures = list(results)
        output_fp = StringIO()

        self.sbol_output.handle_assembly = Mock()

        self.sbol_output.output(results, output_fp)

        self.sbol_output.handle_assembly.assert_has_calls([
            call(sbol_doc_mock, assembly_fixtures[0], 1),
            call(sbol_doc_mock, assembly_fixtures[1], 2)
        ])

        self.assertEqual(output_fp.getvalue(), 'THE SBOL DOC!')

    def test_sanitize_identifier(self):
        """Make sure that sanitize method removes bad characters."""
        bad_examples = [
            ('poop[~1:3]', 'poop__1_3_'),
            ('turdeeee', 'turdeeee'),
            ('pGAL3[~500:2500]', 'pGAL3__500_2500_'),
        ]

        for example_input, expected_output in bad_examples:
            output = self.sbol_output.sanitize_identifier(example_input)
            self.assertEqual(
                output,
                expected_output,
                '"%s" did not result in expected output ""%s". Got "%s"' % (
                    example_input,
                    expected_output,
                    output
                ))

    def test_sbol_handle_assembly(self):
        """Test that we correctly add an assembly to a SBOL document"""
        assembly = self.fixtures.mk_assembly()

        sbol_doc = Document()
        self.sbol_output.handle_assembly(sbol_doc, assembly, 22)

        expected_doc_items = [
            'http://sbols.org/SuperGSL/ComponentDefinition/asm1/1',
            'http://sbols.org/SuperGSL/ComponentDefinition/part-000/1',
            'http://sbols.org/SuperGSL/Sequence/part-000/1',
            'http://sbols.org/SuperGSL/ComponentDefinition/part-001/1',
            'http://sbols.org/SuperGSL/Sequence/part-001/1',
            'http://sbols.org/SuperGSL/Sequence/asm1/1'
        ]
        for expected_uri in expected_doc_items:
            self.assertTrue(
                expected_uri in sbol_doc.SBOLObjects,
                '"%s" not in sbol doc' % expected_uri)
