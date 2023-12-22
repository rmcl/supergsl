from typing import List, Optional
from Bio.Seq import Seq
from supergsl.core.function import SuperGSLFunction
from supergsl.core.sequence import SequenceEntry
from supergsl.core.types.builtin import (
    Collection,
    AminoAcidSequence
)
from supergsl.core.types.protein import Protein
from supergsl.core.types.part import Part
from supergsl.core.types.codon import CodonFrequencyTable

from dnachisel import (
    DnaOptimizationProblem,
    reverse_translate,
    EnforceTranslation,
    AvoidPattern,
    EnforceGCContent,
    MatchTargetCodonUsage
)


class DNAChiselOptimizeFunction(SuperGSLFunction):
    """Run the DNA Chisel Codon Optimization.

    This implementation presently uses a loop where codon-optimized sequences
    are generated sequentially, and each new sequence must differ from all
    previous sequences by at least X nucleotides. This approach was suggested as
    option one in this github issue [3].

    Example sGSL syntax:

        from dnachisel import optimize
        from S288C import codon_frequency

        let aa_sequence = /MRRAGA...RRACC/

        optimize(aa_sequence, codon_frequency, results=5)

    References
        1. https://edinburgh-genome-foundry.github.io/DnaChisel/
        2. https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables/tree/master/codon_usage_data/tables
        3. https://github.com/Edinburgh-Genome-Foundry/DnaChisel/issues/39
    """
    name = 'codon_optimize'
    arguments = [
        ('aa_sequence', (AminoAcidSequence, Protein)),
        #('codon_frequency_table', CodonFrequencyTable),
        ('num_results', int)
    ]
    return_type = Collection

    def execute(self, params : dict):
        """Invoke dnachissel to return matching codon optimized DNA sequence."""
        protein_sequence = params['aa_sequence'].sequence

        protein_name = params['aa_sequence'].identifier

        codon_usage_table = None
        num_results = params['num_results']

        naive_target_sequence = reverse_translate(protein_sequence)
        proposed_sequences : List[SequenceEntry] = []
        for i in range(num_results):
            print('Optimization run %s of %s' % (i, num_results))
            new_sequence = self.create_new_sequence(
                naive_target_sequence=naive_target_sequence,
                codon_usage_table=codon_usage_table,
                existing_sequences=proposed_sequences)
            proposed_sequences.append(new_sequence)

        parts = []
        for idx, sequence_entry in enumerate(proposed_sequences):
            new_part = Part(
                identifier=f'{protein_name}_codon_optimized_{idx}',
                provider=self,
                description=f'Codon optimized version of {protein_name}',
                sequence_entry=sequence_entry)
            parts.append(new_part)

        return Collection(parts)

    def create_new_sequence(
        self,
        naive_target_sequence : str,
        codon_usage_table : Optional[str],
        existing_sequences : List[SequenceEntry]
    ) -> SequenceEntry:
        """Run DNAChisel to create a new codon optimized DNA sequence

        """
        constraints=[
            EnforceTranslation(),
            #EnforceGCContent(mini=0.4, maxi=0.6, window=60),
        ]

        constraints.extend([
            AvoidPattern(str(entry.sequence))
            for entry in existing_sequences
        ])

        problem = DnaOptimizationProblem(
            sequence=naive_target_sequence,
            constraints=constraints,
            objectives=[MatchTargetCodonUsage(species="s_cerevisiae")],
        )

        #print("\nBefore optimization:\n")
        #print(problem.constraints_text_summary())
        #print(problem.objectives_text_summary())

        problem.resolve_constraints(final_check=True)
        problem.optimize()

        #print("\nAfter optimization:\n")
        #print(problem.constraints_text_summary())
        #print(problem.objectives_text_summary())

        sequence_entry = self.sequence_store.add_from_reference(Seq(problem.sequence))
        return sequence_entry
