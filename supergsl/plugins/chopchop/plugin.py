from supergsl.core.plugin import SuperGSLFunction


class ChopChopFunction(SuperGSLFunction):
    """Run the ChopChop CLI tool.

    http://chopchop.cbu.uib.no/
    https://bitbucket.org/valenlab/chopchop
    """

    import_path = 'chopchop'
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

    def get_return_value(self):
        return sgsl_types.PART_LIST

    def execute(self, sgsl_args):
        """
        Generate gRNA for CRISPR/Cas9

        Example sGSL syntax:
            from example import cut
            cut(HO, CAS9, S288C, results=5)
        """
        print('CUT IT UP!')
