import json
import mock
from unittest import TestCase
from supergsl.backend.parts import (
    SliceAndBuildPartSequencePass,
    DNASlice,
    SequencePosition
)
from supergsl.plugins.file_part import FeatureTableWithFastaPartProvider


class PartTestCase(TestCase):

    def setup_genome_provider(self):
        examples = json.load(open('fixtures/s288c_truncated/genome.json'))
        self.examples = {
            ex['description']: ex
            for ex in examples
        }

        settings = {
            'fasta_file_path': 'fixtures/s288c_truncated/genome.fa.gz',
            'feature_file_path': 'fixtures/s288c_truncated/genome.tsv'
        }

        self.part_provider = FeatureTableWithFastaPartProvider('test_genome', settings)

    def test_load_genome_parts(self):
        self.setup_genome_provider()

        part = self.part_provider.get_part('GAL1')
        self.assertEquals(
            part.sequence.find(self.examples['gGAL1']['dna']), 0)

        part = self.part_provider.get_part('HUG1')
        self.assertEquals(
            part.sequence.find(self.examples['gHUG1']['dna']), 0)

    """
    def test_slice_genome_parts(self):
        self.setup_genome_provider()

        for ex_key, ex_data in self.examples.items():
            part = self.part_provider.get_part(ex_key)

            self.assertEquals(
                part.sequence.find(ex_data['dna']), 0, 'Part %s sequence does not match.' % ex_key)
    """
