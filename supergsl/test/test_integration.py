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

    def test_part_slice_gene(self):

        result = self.run_supergsl('''
            from tab.S288C import HO
            gHO''')

        parts = result.get_parts()
        self.assertEquals(len(parts), 1)

        self.assertEquals(
            parts[0].sequence,
            self.expected_sequences.get('gHO').seq)

    def test_part_slice_promoter(self):

        result = self.run_supergsl('''
            from tab.S288C import HO
            pHO''')

        parts = result.get_parts()
        self.assertEquals(len(parts), 1)

        print(parts[0].sequence)
        print()
        print(self.expected_sequences.get('pHO').seq)

        self.assertEquals(
            parts[0].sequence,
            self.expected_sequences.get('pHO').seq)
