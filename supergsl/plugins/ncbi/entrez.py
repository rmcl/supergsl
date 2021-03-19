from .genbank_provider import GenBankFilePartProvider


class EntrezPartProvider(GenBankFilePartProvider):
    """Retrieve genomes and sequences from Entrez-linked public databases.

    Provider Arguments for `supergsl-config.json`:
    `name`: The name of the provider. This is used in your superGSL import statement.
    `provider_class`: For this provider, should be: "supergsl.plugins.entrez.EntrezPartProvider".
    `local_genome_path` (optional): A path on the local filesystem where a local copy
        of the retrieved sequence will be cached.

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
        self.efetch_args = settings['efetch_args']
        self.entrez_email = settings['entrez_email']
        self.genome_path = settings.get(
            'local_genome_path',
            './sgsl-lib/entrez')
        self.allowed_feature_types = set(settings.get('feature_types', [
            'gene'
        ]))

    def get_local_file_path(self):
        os.makedirs(self.genome_path, exist_ok=True)

        # TODO: this makes some strong assumption about the efetch args.
        # refactor to something better later.
        gb_file_path = '%s/%s-%s.gb' % (
            self.genome_path,
            self.efetch_args['db'],
            self.efetch_args['id'],
        )

        return gb_file_path

    def load(self):
        self.retrieve_sequence_file()
        self.genbank_file_path = self.get_local_file_path()
        return super().load()

    def retrieve_sequence_file(self):
        kwargs = self.efetch_args.copy()
        kwargs['rettype'] = 'gb'
        kwargs['retmode'] = 'binary'

        Entrez.email = self.entrez_email
        handle = Entrez.efetch(**kwargs)

        # save to local cache directory
        with open(self.get_local_file_path(), 'w+') as fp:
            fp.write(handle.read())

        handle.close()
