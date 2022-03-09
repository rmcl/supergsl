import gzip
from typing import List, Dict, Tuple, Optional
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from supergsl.core.constants import (
    THREE_PRIME,
    FIVE_PRIME,
    STRAND_WATSON,
    STRAND_CRICK
)
from supergsl.core.exception import PartNotFoundError
from supergsl.core.sequence import SequenceEntry, SequenceAnnotation
from supergsl.core.parts import PartProvider
from supergsl.core.provider import ProviderConfig
from supergsl.core.types.position import Slice, Position
from supergsl.core.types.part import Part


class GenBankFilePartProvider(PartProvider):
    """Load Parts from a Genbank formatted file.

    Provider Arguments for `supergsl-config.json`:
    `name`: The name of the provider. This is used in your superGSL import statement.
    `provider_class`: For this provider, should be: "supergsl.plugins.ncbi.GenBankFilePartProvider".
    `sequence_file_path`: The path on the local filesystem to the genbank file.
    `features` (optional): A list of features to be provided by this provider.
        Default: ['gene', 'CDS']
    `qualifiers` (optional): A list of qualifiers to include as importable parts
        by this provider. Default: ['gene', 'locus_tag', 'product']

    Example of retrieving Adeno-associated virus 4, complete genome (U89790.1)
    {
        "name": "genbank.AAV4",
        "provider_class": "supergsl.plugins.ncbi.GenBankFilePartProvider",
        "sequence_file_path": "/mnt/genomes/nucleotide-U89790.gb.gz"
    }

    More details about GenBank feature tables here: https://www.insdc.org/files/feature_table.html

    If you would like more flexibility in the features that this provider pulls
    out of your GenBane file then you shoould consider subclassing and overriding
    `get_identifier_for_feature`. See method docstring for more information.
    """

    def __init__(self, name, config : ProviderConfig):
        self._provider_name = name

        self.sequence_store = config.sequence_store
        self.genbank_file_path = config.settings['sequence_file_path']
        self.annotation_by_identifier : Dict[str, Tuple[SeqFeature, SequenceEntry]] = {}
        self.loaded = False

        self.file_format = 'genbank'
        self.desired_feature_types = config.settings.get(
            'features',
            ['gene', 'CDS', 'misc_feature'])

        self.desired_feature_qualifiers = config.settings.get(
            'qualifiers',
            ['gene', 'locus_tag', 'product'])

        self.default_part_identifier = config.settings.get(
            'default_part_identifier',
            None)


    def open_gb_file(self, path):
        """Open a genbank file whether gunziped or uncompressed."""
        if path[-2:] == 'gz':
            return gzip.open(path, "rt")

        return open(path, "rt")

    def get_identifier_for_feature(self, feature) -> List[str]:
        """Retrieve part details from a genbank feature.

        This method makes some assumptions about the structure and type of features
        that should be treated as parts and their metadata. You may consider,
        subclassing and overriding this method to handle a your desired genbank
        file.
        """
        identifiers = []
        if feature.type in self.desired_feature_types:
            for qualifier_key in self.desired_feature_qualifiers:
                if qualifier_key in feature.qualifiers:
                    identifiers.extend(feature.qualifiers[qualifier_key])

        return identifiers

    def get_part_from_annotation(
        self,
        identifier : str,
        sequence_annotation : Optional[SequenceAnnotation],
        parent_sequence_entry : SequenceEntry
    ) -> Part:
        """Create a `Part` from the provided feature and reference sequence."""

        if sequence_annotation:
            new_sequence_entry = self.sequence_store.slice_from_annotation(
                parent_sequence_entry,
                sequence_annotation)
            alternative_names = [identifier]

        else:
            new_sequence_entry = parent_sequence_entry
            alternative_names = []

        return Part(
            identifier,
            new_sequence_entry,
            provider=self,
            description=None,
            alternative_names=alternative_names)

    def get_slice_from_feature_location(self, location) -> Slice:
        """Build a slice from a `Bio.SeqFeature.FeatureLocation` object."""
        strand = STRAND_WATSON if location.strand == 1 else STRAND_CRICK
        if strand == STRAND_WATSON:
            return Slice.from_five_prime_indexes(location.start, location.end)

        return Slice(
            Position(location.end * -1, THREE_PRIME, False),
            Position(location.start * -1, THREE_PRIME, False),
            strand=strand)

    def load(self):
        """Open and load features from a genbank file."""
        self.loaded = True
        with self.open_gb_file(self.genbank_file_path) as handle:
            records = SeqIO.parse(handle, self.file_format)

            for record in records:
                sequence_entry = self.sequence_store.add_from_reference(record)

                # Add a record for the entire record
                self.annotation_by_identifier[record.name] = (None, sequence_entry)

                if not self.default_part_identifier:
                    self.default_part_identifier = record.name

                for feature in record.features:
                    annotation_roles = [feature.type]
                    annotation_payload = dict(feature.qualifiers)
                    annotation_payload['type'] = feature.type
                    feature_slice = self.get_slice_from_feature_location(feature.location)

                    new_annotation = SequenceAnnotation(
                        feature_slice,
                        annotation_roles,
                        annotation_payload)
                    sequence_entry.add_annotation(new_annotation)

                    identifiers = self.get_identifier_for_feature(feature)
                    for identifier in identifiers:
                        self.annotation_by_identifier[identifier] = (new_annotation, sequence_entry)

    def list_parts(self) -> List[Part]:
        if not self.loaded:
            self.load()

        return [
            self.get_part(identifier)
            for identifier in self.annotation_by_identifier
        ]

    def get_default_part(self) -> Part:
        """Return the default part for this provider."""
        return self.get_part(self.default_part_identifier)

    def get_part(self, identifier) -> Part:
        """Retrieve a part by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """

        if not self.loaded:
            self.load()

        try:
            sequence_annotation, sequence_entry = self.annotation_by_identifier[identifier]
        except KeyError as error:
            raise PartNotFoundError(
                f'Part not found "{identifier}" in {self.provider_name}.') from error

        return self.get_part_from_annotation(
            identifier,
            sequence_annotation,
            sequence_entry)
