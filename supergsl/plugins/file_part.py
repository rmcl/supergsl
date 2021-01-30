import csv
import gzip
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO

from supergsl.core.constants import THREE_PRIME
from supergsl.core.parts.position import SeqPosition
from supergsl.core.exception import PartLocatorException, PartNotFoundException
from supergsl.core.parts import PartProvider, Part
from supergsl.core.parts.prefix_part import PrefixedSlicePartProviderMixin


class FeatureTableWithFastaPartProvider(PrefixedSlicePartProviderMixin, PartProvider):
    """Access parts provided by fGSL reference genome files.

    The gsl paper [wilson 2008] describes this format as the "Saccharomyces
    Genome Database reference file format". The paper refers to this example:
    http://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab

    These files define the ORF region of genes. GSL assumes that the 500bp preceding
    the ORF is the promoter region and 500 bp downstream are the terminator.

    """

    def __init__(self, name, settings):
        self.name = name
        self.fasta_file_path = settings['fasta_file_path']
        self.feature_file_path = settings['feature_file_path']

    def open_feature_file(self):
        if self.feature_file_path[-2:] == 'gz':
            return gzip.open(self.feature_file_path, "rt")
        else:
            return open(self.feature_file_path, "rt")

    def load(self):
        encoding = guess_type(self.feature_file_path)[1]
        _open_feature_file = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        with _open_feature_file(self.feature_file_path) as handle_fp:
            reader = csv.DictReader(handle_fp, fieldnames=None, delimiter='\t')
            features = list(reader)

            self._genes = {
                feature['gene']: feature
                for feature in features
            }

        encoding = guess_type(self.fasta_file_path)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        with _open(self.fasta_file_path) as fp:
            chromosomes = SeqIO.parse(fp, 'fasta')

            self._sequence_by_chromosome = {
                chromosome.name: chromosome
                for chromosome in chromosomes
            }

    def list_parts(self):
        if not hasattr(self, '_sequence_by_chromosome'):
            self.load()

        return [
            self.get_part(gene_name)
            for gene_name in self._genes.keys()
        ]

    def get_gene(self, gene_name):

        if not hasattr(self, '_sequence_by_chromosome'):
            self.load()

        try:
            reference_feature = self._genes[gene_name]
        except KeyError:
            raise PartNotFoundException('Part not found "%s" in %s.' % (
                gene_name, self.get_provider_name()))

        # Make a copy of the reference feature and modify it to conform
        # to possibly complemented reference sequence
        new_gene_feature = reference_feature.copy()

        chromosome_num = new_gene_feature['chrom#']
        chromosome_sequence = self._sequence_by_chromosome[chromosome_num]

        strand = new_gene_feature['strand']
        if strand == 'C':
            # Translate the position present in the file to be relative to the
            # other strand of DNA (reverse complement).
            reference_sequence = chromosome_sequence.reverse_complement().seq
            reference_len = len(reference_sequence)

            tmp_from = int(new_gene_feature['from'])
            new_gene_feature['from'] = reference_len - int(new_gene_feature['to']) - 1
            new_gene_feature['to'] = reference_len - tmp_from

        else:
            new_gene_feature['from'] = int(new_gene_feature['from'])
            new_gene_feature['to'] = int(new_gene_feature['to']) + 1
            reference_sequence = chromosome_sequence.seq

        return reference_sequence, new_gene_feature

    def _get_alternative_names_from_feature(self, feature):
        alternative_names = set([
            feature['systematic']
        ])
        aliases = feature['aliases'].split(',')
        if len(aliases) > 0 and aliases[0] != '':
            alternative_names.update(aliases)

        return list(alternative_names)

    def get_part(self, identifier):
        """Retrieve a part by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """

        reference_sequence, feature = self.get_gene(identifier)
        alternative_names = self._get_alternative_names_from_feature(feature)

        start = SeqPosition.from_reference(
            x=feature['from'],
            rel_to=THREE_PRIME,
            approximate=False,
            reference=reference_sequence
        )

        print('GET_PART', feature['from'], feature['to'])

        end = start.get_relative_position(
            x=feature['to']-feature['from'])

        part = Part(
            identifier,
            start,
            end,
            provider=self,
            description=feature['Notes'],
            alternative_names=alternative_names)

        return part

    def get_child_part_by_slice(
        self,
        parent_part : Part,
        identifier : str,
        start : SeqPosition,
        end : SeqPosition
    ) -> Part:
        """Return a new part which is the child of the supplied parent."""

        return Part(
            identifier,
            start,
            end,
            provider=self,
            parent_part=parent_part
        )
