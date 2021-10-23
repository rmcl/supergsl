"""Tests for the Assembly module."""
import sys
import unittest
from unittest.mock import Mock, patch

from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.assembly import AssemblyResultOutputFunction

class AssemblyResultOutputFunctionTestCase(unittest.TestCase):
    """Testcases to evaluate the EvaluatePass class."""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.function_config = self.fixtures.mk_function_config_object()

    def test_open_output_fp_stdout(self):
        """Test that open_output_fp returns a filepointer to stdout."""
        output = AssemblyResultOutputFunction(self.function_config)

        with output.open_output_fp('-') as fp:
            self.assertEqual(fp, sys.stdout)

    @patch('supergsl.core.assembly.open')
    def test_open_output_fp_file_open(self, open_mock):
        """Test that open_output_fp opens a file and then closes it when done."""
        open_mock.return_value = file_handle_mock = Mock()
        output = AssemblyResultOutputFunction(self.function_config)

        with output.open_output_fp('test.txt') as output_fp:
            output_fp.write('hello world!')

        file_handle_mock.write.assert_called_once_with('hello world!')
        file_handle_mock.close.assert_called_once_with()

    def test_execute_calls_expected_methods(self):
        """The execute method should call output with correct arguments."""
        output = AssemblyResultOutputFunction(self.function_config)
        output.output = Mock()
        assemblies = ['hello', 'assemblies']
        output.execute({
            'assemblies': assemblies,
            'filename': '-'
        })

        output.output.assert_called_once_with(
            assemblies,
            sys.stdout
        )
