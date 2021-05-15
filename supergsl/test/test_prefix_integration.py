from unittest import TestCase
from Bio import SeqIO
from supergsl.core.pipeline import CompilerPipeline
from supergsl.test.fixtures import SuperGSLIntegrationFixtures
from supergsl.test.fixtures.utils import TestOutputAstPass

class SuperGSLIntegrationTestCases(TestCase):

    def setUp(self):
        self.expected_sequences = SeqIO.index(
            'supergsl/test/expected_sequences.fasta', 'fasta')

        self.fixtures = SuperGSLIntegrationFixtures()
        self.compiler_settings = self.fixtures.get_supergsl_settings()

    def run_supergsl(self, source_code):
        pipeline = CompilerPipeline(self.compiler_settings)
        pipeline.compile(source_code)

        return pipeline

    def test_part_prefix_slice_locus(self):

        gsl_template = '''
            from truncated.S288C import HO
            let test_part = %s'''

        # Currently ommitting dHO and tHO. I think that the primer designer may
        # shifting these in fGSL to make better parts so will revisit once we get
        # to implementing that
        things_to_test = [
            'gHO',
            'pHO',
            'uHO',
        ]

        for part_name in things_to_test:
            result = self.run_supergsl(gsl_template % part_name)

            symbol_table = result.get_symbol_table()
            part = symbol_table.lookup('test_part')
            self.assertEqual(
                part.get_sequence().seq,
                self.expected_sequences.get(part_name).seq,
                '%s sequence does not match expection' % part_name)
