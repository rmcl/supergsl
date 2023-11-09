import csv
import gzip
from typing import Dict, List, Tuple, Any, TextIO
from mimetypes import guess_type
from Bio import SeqIO
from Bio.Seq import Seq

from supergsl.core.constants import (
    FIVE_PRIME,
    THREE_PRIME,
    STRAND_WATSON,
    STRAND_CRICK
)
from supergsl.core.types.position import Slice, Position
from supergsl.core.sequence import SequenceEntry
from supergsl.core.exception import PartNotFoundError
from supergsl.core.types.part import Part
from supergsl.core.parts import PartProvider
from supergsl.core.provider import ProviderConfig
from supergsl.core.parts.prefix_part import PrefixedSlicePartProviderMixin
from supergsl.plugins.pydna.primers import ExtractionPrimerBuilder

class FeatureTableWithFastaPartProvider(PrefixedSlicePartProviderMixin, PartProvider):
    """Access parts provided by fGSL reference genome files.

    The gsl paper [wilson 2008] describes this format as the "Saccharomyces
    Genome Database reference file format". The paper refers to this example:
    http://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab

    These files define the ORF region of genes. GSL assumes that the 500bp preceding
    the ORF is the promoter region and 500 bp downstream are the terminator.

    """

    def __init__(self, name : str, config : ProviderConfig):
        self._provider_name = name
        self.sequence_store = config.sequence_store

        settings = config.settings
        self.fasta_file_path : str = settings['fasta_file_path']
        self.feature_file_path : str = settings['feature_file_path']
        self._cached_parts : Dict[str, Part] = {}
        self.primer_builder = ExtractionPrimerBuilder(self.sequence_store)

    def load(self) -> None:
        encoding = guess_type(self.feature_file_path)[1]

        def _open(file_path : str) -> TextIO:
            encoding = guess_type(file_path)[1]
            if encoding == 'gzip':
                return gzip.open(file_path, mode='rt')
            else:
                return open(file_path, mode='rt')

        with _open(self.feature_file_path) as handle_fp:
            reader = csv.DictReader(handle_fp, fieldnames=None, delimiter='\t')
            features = list(reader)

            self._genes = {
                feature['gene']: feature
                for feature in features
            }

        with _open(self.fasta_file_path) as fp:
            chromosomes = SeqIO.parse(fp, 'fasta')

            self._sequence_by_chromosome = {}
            for chromosome in chromosomes:

                chromosome_entry = self.sequence_store.add_from_reference(chromosome.seq)
                # Todo: Add annotations for this chromosome_entry

                self._sequence_by_chromosome[chromosome.name] = chromosome_entry

    def list_parts(self):
        if not hasattr(self, '_sequence_by_chromosome'):
            self.load()

        return [
            self.get_part(gene_name)
            for gene_name in self._genes.keys()
        ]

    def get_gene(self, gene_name : str) -> Tuple[SequenceEntry, dict]:

        if not hasattr(self, '_sequence_by_chromosome'):
            self.load()

        try:
            reference_feature = self._genes[gene_name]
        except KeyError:
            raise PartNotFoundError('Part not found "%s" in %s.' % (
                gene_name, self.provider_name))

        # Make a copy of the reference feature and modify it to conform
        # to possibly complemented reference sequence
        new_gene_feature : Dict[str, Any] = reference_feature.copy()

        chromosome_num = new_gene_feature['chrom#']
        chromosome_sequence_entry = self._sequence_by_chromosome[chromosome_num]

        strand = new_gene_feature['strand']
        if strand == 'C':
            # This gene is on the reverse strand.
            sequence_len = len(chromosome_sequence_entry.sequence)
            new_gene_feature['from'] = sequence_len - int(reference_feature['to']) - 1
            new_gene_feature['to'] = sequence_len - int(reference_feature['from'])
            strand = STRAND_CRICK

        else:
            # This gene is on the forward strand.
            new_gene_feature['from'] = int(reference_feature['from'])
            new_gene_feature['to'] = int(reference_feature['to']) + 1
            strand = STRAND_WATSON


        part_slice = Slice(
            Position(new_gene_feature['from'], FIVE_PRIME, False),
            Position(new_gene_feature['to'], FIVE_PRIME, False),
            strand=strand)

        new_entry = self.sequence_store.slice(
            chromosome_sequence_entry,
            part_slice)

        return new_entry, new_gene_feature

    def _get_alternative_names_from_feature(self, feature : dict) -> List[str]:
        alternative_names = set([
            feature['systematic']
        ])
        aliases = feature['aliases'].split(',')
        if len(aliases) > 0 and aliases[0] != '':
            alternative_names.update(aliases)

        return list(alternative_names)

    def get_part(self, identifier : str) -> Part:
        """Retrieve a part by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """
        try:
            return self._cached_parts[identifier]
        except KeyError:
            pass

        sequence_entry, feature = self.get_gene(identifier)
        alternative_names = self._get_alternative_names_from_feature(feature)


        part = Part(
            identifier,
            sequence_entry,
            provider=self,
            description=feature['Notes'],
            alternative_names=alternative_names)

        # Build primers for this part.
        #self.primer_builder.build_primers_for_part(part)

        self._cached_parts[identifier] = part
        return part

    def get_child_part_by_slice(
        self,
        parent_part : Part,
        identifier : str,
        part_slice : Slice,
    ) -> Part:
        """Return a new part which is the child of the supplied parent."""

        new_sequence_entry = self.sequence_store.slice(
            parent_part.sequence_entry,
            part_slice)

        return Part(
            identifier,
            sequence_entry=new_sequence_entry,
            provider=self,
            parent_part=parent_part
        )
