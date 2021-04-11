from supergsl.core.function import SuperGSLFunction
from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.ast import Assembly
from supergsl.core.types import NucleotideSequence



class ChopChopFunction(SuperGSLPlugin):
    """Run the ChopChop CLI tool.

    http://chopchop.cbu.uib.no/
    https://bitbucket.org/valenlab/chopchop
    """

    name = 'cut'

    def get_help(self):
        return self.execute.__docstr__

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

    def register(self, symbol_table, compiler_settings):
        """Register functions provide by chopchop."""
        provider_table = symbol_table.nested_scope('imports')
        provider_table.insert('chopchop', ChopChopFunction())
