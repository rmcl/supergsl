import unittest
from Bio.Seq import Seq
from supergsl.core.function import SuperGSLFunction
from supergsl.core.types.builtin import AminoAcidSequence, NucleotideSequence
from supergsl.core.exception import FunctionInvokeError

class SuperGSLFunctionTestCase(unittest.TestCase):
    """Test that SuperGSL functions can process arguments and be invoked."""

    class TestFunction(SuperGSLFunction):
        arguments = [
            ('awesome1', int),
            ('arg2', AminoAcidSequence)
        ]

    def setUp(self):
        compiler_settings = {}
        self.function = SuperGSLFunctionTestCase.TestFunction(compiler_settings)

    def test_build_argument_map_different_num_positional_args(self):
        """Raise exception when the number of positional arguments differs from expected."""

        positional_args = []
        child_arguments = []

        self.assertRaises(
            FunctionInvokeError,
            self.function.build_argument_map,
            positional_args,
            child_arguments)

    def test_build_argument_map_bad_argument_type(self):
        """Raise exception when one of the arguments does not match expected type."""

        positional_args = [
            55,
            NucleotideSequence(Seq('ATGC'))
        ]
        child_arguments = []

        self.assertRaises(
            FunctionInvokeError,
            self.function.build_argument_map,
            positional_args,
            child_arguments)

    def test_build_argument_map(self):
        """Return the correct dictionary of argument names to their passed values."""
        positional_args = [
            55,
            AminoAcidSequence(Seq('MAATT*'))
        ]
        child_arguments = []

        result = self.function.build_argument_map(positional_args, child_arguments)

        self.assertEqual(result, {
            'arg2': positional_args[1],
            'awesome1': 55
        })
