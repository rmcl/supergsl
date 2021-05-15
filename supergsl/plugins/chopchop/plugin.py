from supergsl.core.function import SuperGSLFunction, SuperGSLFunctionDeclaration
from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunction
from supergsl.core.ast import Assembly
from supergsl.core.types.builtin import NucleotideSequence



class ChopChopFunction(SuperGSLFunction):
    """Run the ChopChop CLI tool.

    http://chopchop.cbu.uib.no/
    https://bitbucket.org/valenlab/chopchop
    """

    name = 'cut'

    def get_arguments(self):
        return [
            argument('target_gene', sgsl_types.PART),
            argument('target_genome', sgsl_types.GENOME),
            argument('target_genome', sgsl_types.GENOME),
            argument('num_results', int)
        ]

    def get_return_type(self):
        return NucleotideConstant

    def execute(self, sgsl_args, child_nodes=None):
        """
        Generate gRNA for CRISPR/Cas9

        Example sGSL syntax:
            from example import cut
            cut(HO, CAS9, S288C, results=5)
        """
        self.invoke()

        print('CUT IT UP!')
        return NucleotideSequence('TTA')


class ChopChopPlugin(SuperGSLPlugin):

    def register(self, compiler_settings : dict):
        """Register functions provide by chopchop."""
        self.register_function(
            'chopchop',
            'chopchop',
            SuperGSLFunctionDeclaration(ChopChopFunction, compiler_settings))
