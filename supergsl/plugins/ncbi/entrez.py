import os
from Bio import Entrez
from supergsl.core.exception import ConfigurationError
from .genbank_provider import GenBankFilePartProvider


class EntrezPartProvider(GenBankFilePartProvider):
    """Retrieve genomes and sequences from Entrez-linked public databases.

    Provider Arguments for `supergsl-config.json`:
    `name`: The name of the provider. This is used in your superGSL import statement.
    `provider_class`: For this provider, should be: "supergsl.plugins.entrez.EntrezPartProvider".
    `local_genome_path` (optional): A path on the local filesystem where a local copy
        of the retrieved sequence will be cached.
    `cache_records` (optional): Boolean flag, if False re-retrieve the sequence file
        from Entrez every time. Defaults to True.

    Example of retrieving Adeno-associated virus 4, complete genome (U89790.1)
    {
        "name": "entrez.AAV4",
        "provider_class": "supergsl.plugins.entrez.EntrezPartProvider",
        "efetch_args": {
            "db": "nucleotide",
            "id": "U89790"
        }
    }

    """

    def __init__(self, name, settings):
        self.name = name
        self.efetch_args = settings.get('efetch_args', None)
        self.entrez_email = settings.get('entrez_email', None)
        self.use_cached_records = settings.get('cache_records', True)
        self.genome_path = settings.get(
            'local_genome_path',
            './sgsl-lib/entrez/')
        self.allowed_feature_types = set(settings.get('feature_types', [
            'gene'
        ]))

        if not self.efetch_args or not self.entrez_email:
            raise ConfigurationError(
                'Entrez provider not correctly configured. `efetch_args and` '
                '`entrez_email` must be specified.')

    def get_local_file_path(self):
        """Get the path to the local cache genbank file for this Entrez entry."""
        os.makedirs(self.genome_path, exist_ok=True)

        # TODO: this makes some strong assumption about the efetch args.
        # refactor to something better later.

        filename = '%s-%s.gb' % (
            self.efetch_args['db'],
            self.efetch_args['id'],
        )
        return os.path.join(self.genome_path, filename)

    def cache_record_exists(self):
        """Determine if a cache exists for a retrieved record.

        This method only checks if file exists. Consider better hash checking in
        the future if supported by Entrez.
        """
        return os.path.exists(self.get_local_file_path())

    def load(self):
        """Retrieve record for Entrez and load identifiers to allow their potential import."""
        if not (self.use_cached_records and self.cache_record_exists()):
            self.retrieve_sequence_file()

        self.genbank_file_path = self.get_local_file_path()
        return super().load()

    def _open_local_gb_file(self):
        return open(self.get_local_file_path(), 'w+')

    def retrieve_sequence_file(self):
        """Retrieve a genbank file from Entrez servers."""
        kwargs = self.efetch_args.copy()
        kwargs['rettype'] = 'gb'
        kwargs['retmode'] = 'binary'

        Entrez.email = self.entrez_email
        handle = Entrez.efetch(**kwargs)

        with self._open_local_gb_file() as local_gb_fp:
            local_gb_fp.write(handle.read())

        handle.close()
