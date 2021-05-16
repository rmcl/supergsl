import csv
import gzip
from typing import Dict, List, Tuple, Any, TextIO
from mimetypes import guess_type
from Bio import SeqIO
from Bio.Seq import Seq

from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.position import SeqPosition
from supergsl.core.exception import PartNotFoundError
from supergsl.core.types.part import Part
from supergsl.core.parts import PartProvider
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

    def __init__(self, name : str, settings : dict):
        self.name = name
        self.fasta_file_path : str = settings['fasta_file_path']
        self.feature_file_path : str = settings['feature_file_path']
        self._cached_parts : Dict[str, Part] = {}
        self.primer_builder = ExtractionPrimerBuilder()

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

    def get_gene(self, gene_name : str) -> Tuple[Seq, dict]:

        if not hasattr(self, '_sequence_by_chromosome'):
            self.load()

        try:
            reference_feature = self._genes[gene_name]
        except KeyError:
            raise PartNotFoundError('Part not found "%s" in %s.' % (
                gene_name, self.get_provider_name()))

        # Make a copy of the reference feature and modify it to conform
        # to possibly complemented reference sequence
        new_gene_feature : Dict[str, Any] = reference_feature.copy()

        chromosome_num = new_gene_feature['chrom#']
        chromosome_sequence = self._sequence_by_chromosome[chromosome_num]

        strand = new_gene_feature['strand']
        if strand == 'C':
            # Translate the position present in the file to be relative to the
            # other strand of DNA (reverse complement).
            reference_sequence = chromosome_sequence.reverse_complement().seq
            reference_len = len(reference_sequence)

            new_gene_feature['from'] = reference_len - int(reference_feature['to']) - 1
            new_gene_feature['to'] = reference_len - int(reference_feature['from'])

        else:
            new_gene_feature['from'] = int(reference_feature['from'])
            new_gene_feature['to'] = int(reference_feature['to']) + 1
            reference_sequence = chromosome_sequence.seq

        return reference_sequence, new_gene_feature

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

        reference_sequence, feature = self.get_gene(identifier)
        alternative_names = self._get_alternative_names_from_feature(feature)

        start = SeqPosition.from_reference(
            x=feature['from'],
            rel_to=THREE_PRIME,
            approximate=False,
            reference=reference_sequence
        )

        end = start.get_relative_position(
            x=feature['to']-feature['from'])

        part = Part(
            identifier,
            start,
            end,
            provider=self,
            description=feature['Notes'],
            alternative_names=alternative_names)

        # Build primers for this part.
        self.primer_builder.build_primers_for_part(part)

        self._cached_parts[identifier] = part
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
