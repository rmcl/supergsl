import gzip
from typing import List
from Bio import SeqIO
from supergsl.core.constants import THREE_PRIME
from supergsl.core.exception import PartNotFoundError
from supergsl.core.parts import PartProvider
from supergsl.core.types.part import Part
from supergsl.core.types.position import SeqPosition



class GenBankFilePartProvider(PartProvider):
    """Load Parts from a Genbank formatted file.

    Provider Arguments for `supergsl-config.json`:
    `name`: The name of the provider. This is used in your superGSL import statement.
    `provider_class`: For this provider, should be: "supergsl.plugins.ncbi.GenBankFilePartProvider".
    `sequence_file_path`: The path on the local filesystem to the genbank file.

    Example of retrieving Adeno-associated virus 4, complete genome (U89790.1)
    {
        "name": "genbank.AAV4",
        "provider_class": "supergsl.plugins.ncbi.GenBankFilePartProvider",
        "sequence_file_path": "/mnt/genomes/nucleotide-U89790.gb.gz"
    }

    If you don't like the features that this provider pulls out of your GenBank
    file then you shoould consider subclassing and overriding
    `get_identifier_for_feature`. See method docstring for more information.
    """

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

    def get_identifier_for_feature(self, feature) -> List[str]:
        """Retrieve part details from a genbank feature.

        This method makes some assumptions about the structure and type of features
        that should be treated as parts and their metadata. You may consider,
        subclassing and overriding this method to handle a your desired genbank
        file.
        """
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
            records = SeqIO.parse(handle, 'genbank')

            for record in records:
                for feature in record.features:
                    identifiers = self.get_identifier_for_feature(feature)
                    for identifier in identifiers:
                        self.features_by_identifier[identifier] = (feature, record)

    def list_parts(self) -> List[Part]:
        if not self.loaded:
            self.load()

        return [
            self.get_part(identifier)
            for identifier in self.features_by_identifier.keys()
        ]

    def get_part(self, identifier) -> Part:
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
            raise PartNotFoundError('Part not found "%s" in %s.' % (
                identifier, self.get_provider_name()))

        return self.get_part_from_feature(
            identifier,
            feature,
            parent_record)
