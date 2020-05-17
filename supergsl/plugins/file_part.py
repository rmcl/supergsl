import gzip
from Bio import SeqIO

from supergsl.core.exception import PartLocatorException
from supergsl.backend.parts import PartProvider, Part


class GenbankFilePartProvider(PartProvider):

    def __init__(self, name, settings):
        self.name = name
        self.genbank_file_path = settings['sequence_file_path']

        self.load()

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
        try:
            feature, parent_record = self.features_by_gene_name[identifier]
        except KeyError:
            raise PartLocatorException('Part not found "%s" in %s.' % (identifier, self.get_provider_name()))

        location = feature.location
        sequence = parent_record[location.start:location.end]

        part = Part(identifier, sequence)
        part.feature = feature
        return part
