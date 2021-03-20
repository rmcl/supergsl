import os
import gzip
from Bio import Entrez, SeqIO
from supergsl.core.constants import THREE_PRIME
from supergsl.core.parts import PartProvider, Part, SeqPosition



class GenBankFilePartProvider(PartProvider):
    """Load Parts from a Genbank formatted file."""

    def __init__(self, name, settings):
        self.name = name
        self.genbank_file_path = settings['sequence_file_path']
        self.features_by_identifier = {}
        self.loaded = False

    def open_gb_file(self, path):
        """Open a genbank file whether gunziped or uncompressed."""
        if path[-2:] == 'gz':
            return gzip.open(path, "rt")
        else:
            return open(path, "rt")

    def get_identifier_for_feature(self, feature) -> str:
        identifiers = []
        if feature.type in ['gene', 'CDS']:
            desired_qualifiers = ['gene', 'locus_tag', 'product']
            for qualifier_key in desired_qualifiers:
                if qualifier_key in feature.qualifiers:
                    identifiers.extend(feature.qualifiers[qualifier_key])

        return identifiers

    def get_part_from_feature(self, identifier, feature, parent_sequence_record) -> Part:
        """Create a `Part` from the provided feature and reference sequence."""

        location = feature.location
        start = SeqPosition.from_reference(
            x=location.start,
            rel_to=THREE_PRIME,
            approximate=False,
            reference=parent_sequence_record.seq)

        end = start.get_relative_position(
            x=location.end-location.start)

        alternative_names = self.get_identifier_for_feature(feature)
        return Part(
            identifier,
            start,
            end,
            provider=self,
            description=None,
            alternative_names=alternative_names)


    def load(self):
        """Open and load features from a genbank file."""
        self.loaded = True
        with self.open_gb_file(self.genbank_file_path) as handle:
            records = SeqIO.parse(
                handle,
                'genbank')

            for record in records:
                for feature in record.features:
                    identifiers = self.get_identifier_for_feature(feature)
                    for identifier in identifiers:
                        self.features_by_identifier[identifier] = (feature, record)

    def list_parts(self):
        if not self.loaded:
            self.load()

        return [
            self.get_part(identifier)
            for identifier in self.features_by_identifier.keys()
        ]

    def get_part(self, identifier):
        """Retrieve a part by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """

        if not self.loaded:
            self.load()

        try:
            feature, parent_record = self.features_by_identifier[identifier]
        except KeyError:
            raise PartLocatorException('Part not found "%s" in %s.' % (identifier, self.get_provider_name()))


        part = self.get_part_from_feature(
            identifier,
            feature,
            parent_record)
        return part
