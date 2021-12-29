import unittest
from Bio.Seq import Seq
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.function import SuperGSLFunction
from supergsl.core.types.builtin import AminoAcidSequence, NucleotideSequence
from supergsl.core.exception import FunctionInvokeError

class SuperGSLFunctionTestCase(unittest.TestCase):
    """Test that SuperGSL functions can process arguments and be invoked."""

    class TestFunction(SuperGSLFunction):
        name = 'party'
        arguments = [
            ('awesome1', int),
            ('arg2', AminoAcidSequence)
        ]

        return_type = NucleotideSequence

        def execute(self, params : dict):
            sequence_entry = self.sequence_store.add_from_reference(Seq('ATG'))
            return NucleotideSequence(sequence_entry)

    def setUp(self):
        compiler_settings = {}
        self.fixtures = SuperGSLCoreFixtures()
        self.function_config = self.fixtures.mk_provider_config(compiler_settings)
        self.function = SuperGSLFunctionTestCase.TestFunction(self.function_config)

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

    def test_argument_named_children_disallowed(self):
        """Raise exception if the function defines a argument "children"."""

        class TestFunction2(SuperGSLFunction):
            """Test SuperGSLFunction definition."""
            arguments = [
                ('children', int),
            ]

        function = TestFunction2(self.function_config)

        self.assertRaises(
            FunctionInvokeError,
            function.build_argument_map,
            [25],
            [])


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

    def test_check_function_result_return_value_mismatch(self):
        """Raise exception if the return value doesn't match expected type."""
        result = 55
        self.assertRaises(
            FunctionInvokeError,
            self.function.check_function_result,
            result)

    def test_evaluate_arguments_and_execute(self):
        positional_arguments = [
            55,
            AminoAcidSequence(Seq('MAATT*'))
        ]
        child_arguments = []
        result = self.function.evaluate_arguments_and_execute(
            positional_arguments, child_arguments)

        self.assertEqual(result.sequence, 'ATG')
