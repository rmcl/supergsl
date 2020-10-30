from unittest import TestCase
from Bio import SeqIO
from supergsl.core.pipeline import CompilerPipeline
from supergsl.core.output import TestOutputProvider

class SuperGSLIntegrationTestCases(TestCase):

    def setUp(self):
        self.expected_sequences = SeqIO.index(
            'supergsl/test/expected_sequences.fasta', 'fasta')

    def run_supergsl(self, source_code):
        pipeline = CompilerPipeline()
        ast = pipeline.compile(source_code)

        output = TestOutputProvider()
        output.perform(ast)

        return output

    def test_part_slice_HO_locus(self):

        gsl_template = '''
            from tab.S288C import HO
            %s'''

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

            parts = result.get_parts()
            self.assertEquals(len(parts), 1)
            self.assertEquals(
                parts[0].get_sequence(),
                self.expected_sequences.get(part_name).seq,
                '%s sequence does not match expection' % part_name)
