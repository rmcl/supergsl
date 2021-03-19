import os
import gzip
from Bio import Entrez, SeqIO
from supergsl.core.parts import PartProvider


class GenBankFilePartProvider(PartProvider):
    """Load Parts from a Genbank formatted file."""

    def __init__(self, name, settings):
        self.name = name
        self.genbank_file_path = settings['sequence_file_path']
        self.allowed_feature_types = set(settings.get('feature_types', [
            'gene'
        ]))

    def open_gb_file(self, path):
        """Open a genbank file whether gunziped or uncompressed."""
        if path[-2:] == 'gz':
            return gzip.open(path, "rt")
        else:
            return open(path, "rt")

    def load(self):
        """Open and load features from a genbank file."""
        with self.open_gb_file(self.genbank_file_path) as handle:
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
                        if feature.type in self.allowed_feature_types and 'gene' in feature.qualifiers
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

    def list_parts(self):
        if not hasattr(self, 'features_by_gene_name'):
            self.load()

        return [
            self.get_part(gene_name)
            for gene_name in self.features_by_gene_name.keys()
        ]

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
