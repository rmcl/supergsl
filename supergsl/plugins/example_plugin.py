class ExamplePlugin(SuperGSLPlugin):

    # Types of plugins
    #    - SuperGSL Function


    def register(self, context):
        cut_func = FunctionRegistration(
            'cut', self.cut, return_value=sgsl_types.PART_LIST, help=self.cut.__docstr__)

        cut_fun.add_argument('target_gene', sgsl_types.PART)
        cut_fun.add_argument('target_genome', sgsl_types.GENOME)
        cut_fun.add_argument('target_genome', sgsl_types.GENOME)
        cut_fun.add_argument('num_results', int)

    def cut(self, sgsl_args):
        """
        Generate gRNA for CRISPR/Cas9

        Example sGSL syntax:
            from example import cut
            cut(HO, CAS9, S288C, results=5)
        """
