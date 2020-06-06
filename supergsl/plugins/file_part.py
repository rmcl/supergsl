import csv
import gzip
from Bio import SeqIO

from supergsl.core.exception import PartLocatorException, PartNotFoundException
from supergsl.backend.parts import PartProvider, Part


class FeatureTableWithFastaPartProvider(PartProvider):
    """This provider matches the Reference Genome files that GSL 1.0 utilizes."""

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
        with self.open_feature_file() as handle_fp:
            reader = csv.DictReader(handle_fp, fieldnames=None, delimiter='\t')
            features = list(reader)

            self._genes = {
                feature['gene']: feature
                for feature in features
            }

        self._sequence_by_chromosome = list(SeqIO.parse(self.fasta_file_path, 'fasta'))

    def get_gene(self, gene_name):

        if not hasattr(self, '_sequence_by_chromosome'):
            self.load()

        try:
            feature = self._genes[gene_name]
        except KeyError:
            raise PartNotFoundException('Part not found "%s" in %s.' % (
                gene_name, self.get_provider_name()))

        chromosome_num = int(feature['chrom#'])
        chromosome_sequence = self._sequence_by_chromosome[chromosome_num - 1]

        loc = (
            int(feature['from']),
            int(feature['to']) + 1 # GSL uses non-zero relative indexes!! Do we want to conform to this insanity???
        )

        print(gene_name, 'chrom', chromosome_sequence.id, loc, loc[1] - loc[0])

        strand = feature['strand']
        if strand == 'C':
            seq = chromosome_sequence[loc[0]:loc[1]].reverse_complement().seq
        else:
            seq = chromosome_sequence[loc[0]:loc[1]].seq

        return seq, feature

    def get_part(self, identifier):
        """Retrieve a part by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """
        sequence, feature = self.get_gene(identifier)
        part = Part(identifier, sequence)
        part.feature = feature
        return part


class GenbankFilePartProvider(PartProvider):

    def __init__(self, name, settings):
        self.name = name
        self.genbank_file_path = settings['sequence_file_path']

    def load(self):
        with gzip.open(self.genbank_file_path, "rt") as handle:
            records = SeqIO.parse(
                handle,
                'genbank'
            )

            features_by_gene_name = []
            features_by_systematic_name = []
            for record in records:
                features_by_gene_name += [
                    (feature.qualifiers['gene'][0], (feature, record))
                    for feature in record.features
                    if feature.type == 'gene' and 'gene' in feature.qualifiers
                ]

            """
                features_by_systematic_name += [
                    (feature.qualifiers['locus_tag'][0], feature)
                    for feature in record.features
                    if feature.type == 'gene'
                ]
            """

            self.features_by_gene_name = dict(features_by_gene_name)
            #self.features_by_systematic_name = dict(features_by_systematic_name)

    def get_part(self, identifier):
        """Retrieve a part by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """

        if not hasattr(self, 'features_by_gene_name'):
            self.load()

        try:
            feature, parent_record = self.features_by_gene_name[identifier]
        except KeyError:
            raise PartLocatorException('Part not found "%s" in %s.' % (identifier, self.get_provider_name()))

        location = feature.location
        sequence = parent_record[location.start:location.end]

        part = Part(identifier, sequence)
        part.feature = feature
        return part
