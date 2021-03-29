from supergsl.core.types import SuperGSLEnum
from supergsl.core.function import SuperGSLFunction
from supergsl.core.types import NucleotideSequence, AminoAcidSequence, CodonFrequencyTable

from dnachisel import (
    DnaOptimizationProblem,
    reverse_translate,
    EnforceTranslation,
    AvoidPattern,
    EnforceGCContent,
    MatchTargetCodonUsage
)


class DNAChiselFunction(SuperGSLFunction):
    """Run the DNA Chisel Codon Optimization.

    This implementation presently uses a loop where codon-optimized sequences
    are generated sequentially, and each new sequence must differ from all
    previous sequences by at least X nucleotides. This approach was suggested as
    option one in this github issue [3].

    Example sGSL syntax:
        from dnachisel import optimize
        from s288c import codon_frequency

        let aa_sequence = /MRRAGA...RRACC/

        optimize(aa_sequence, codon_frequency, results=5)

    References
        * https://edinburgh-genome-foundry.github.io/DnaChisel/
        * https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables/tree/master/codon_usage_data/tables
        * https://github.com/Edinburgh-Genome-Foundry/DnaChisel/issues/39
    """

    import_path = 'dnachisel'
    name = 'codon_optimize'
    arguments = [
        ('aa_sequence', AminoAcidSequence),
        ('codon_frequency_table', CodonFrequencyTable),
        ('num_results', int)
    ]
    return_type = NucleotideSequence

    def execute(self, sgsl_args, child_nodes=None):
        """Invoke dnachissel to return matching codon optimized DNA sequence."""

        protein = sgsl_args['aa_sequence']
        naive_target_sequence = reverse_translate(protein)
        codon_usage_table = None
        num_results = sgsl_args.num_results

        proposed_sequences = []
        for i in range(num_results):
            print('Optimization run %s of %s' % (i, num_results))
            new_sequence = self.create_new_sequence(
                naive_target_sequence=naive_target_sequence,
                codon_usage_table=codon_usage_table,
                existing_sequences=proposed_sequences)
            proposed_sequences.append(new_sequence)

        return [
            NucleotideSequence(sequence)
            for sequence in proposed_sequences
        ]

    def create_new_sequence(
        self,
        naive_target_sequence,
        codon_usage_table,
        existing_sequences
    ):
        """Run DNAChisel to create a new codon optimized DNA sequence

        """
        constraints=[
            EnforceTranslation(),
            #EnforceGCContent(mini=0.4, maxi=0.6, window=60),
        ]

        constraints.extend([
            AvoidPattern(sequence)
            for sequence in existing_sequences
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

        return problem.sequence
