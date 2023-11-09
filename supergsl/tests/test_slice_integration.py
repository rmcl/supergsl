from unittest import TestCase
from Bio import SeqIO
from supergsl.lang.pipeline import CompilerPipeline
from supergsl.tests.fixtures import SuperGSLIntegrationFixtures
from supergsl.tests.fixtures.utils import TestOutputAstPass

from supergsl.utils import display_symbol_table


class SuperGSLIntegrationTestCases(TestCase):

    def setUp(self):
        self.maxDiff = None
        self.expected_sequences = SeqIO.index(
            'supergsl/tests/expected_sequences.fasta', 'fasta')

        self.fixtures = SuperGSLIntegrationFixtures()
        self.compiler_settings = self.fixtures.get_supergsl_settings()

    def run_supergsl(self, source_code):
        pipeline = CompilerPipeline(self.compiler_settings)
        ast = pipeline.compile(source_code)

        return pipeline

    def test_part_slice_notation(self):

        gsl_template = '''
            from truncated.S288C import HMG1
            let test_part = %s'''

        things_to_test = [
            'gHMG1[0:100S]',
            'gHMG1[-500E:100E]',
            'gHMG1[2087S:200E]',
            'gHMG1[1586:200E]'
        ]

        for part_name in things_to_test:
            result = self.run_supergsl(gsl_template % part_name)

            part = result.symbols['test_part']
            expected = self.expected_sequences.get(part_name, None)
            assert expected is not None, 'Could not find sequence %s' % (part_name)

            self.assertEqual(
                part.sequence,
                self.expected_sequences.get(part_name).seq,
                '%s sequence does not match expection' % part_name)
