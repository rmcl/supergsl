from supergsl.plugins.docker import DockerFunction
from supergsl.core.types.part import Part
from supergsl.core.types.assembly import Assembly
from supergsl.core.types.builtin import NucleotideSequence



class GuideRNABuilderFunction(DockerFunction):
    """Run the ChopChop CLI tool to create a cas9 gRNA for cutting a specific locus.

    http://chopchop.cbu.uib.no/
    https://bitbucket.org/valenlab/chopchop
    """

    arguments = [
        ('target_gene', Part),
        #('target_genome', Genome),
        ('num_results', int)
    ]
    return_type = NucleotideSequence
    image_tag = 'chopchop:latest'


    def get_docker_command(self, params):
        return 'hello'
    """
    def execute(self, params : dict):
        ""
        Generate gRNA for CRISPR/Cas9

        Example sGSL syntax:
            from example import cut
            cut(HO, CAS9, S288C, results=5)
        ""
        self.invoke()

        print('CUT IT UP!')
        return NucleotideSequence('TTA')
    """
