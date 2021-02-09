from Bio.Seq import Seq
from supergsl.core.ast import Assembly, Part as AstPart
from supergsl.core.constants import THREE_PRIME
from supergsl.core.parts import Part, SeqPosition

class SuperGSLIntegrationFixtures(object):

    def get_supergsl_settings(self):
        return {
            "part_providers": [
                {
                    "name": "truncated.S288C",
                    "provider_class": "supergsl.plugins.file_part.FeatureTableWithFastaPartProvider",
                    "fasta_file_path": "supergsl/test/fixtures/S288C_truncated/genome.fa.gz",
                    "feature_file_path": "supergsl/test/fixtures/S288C_truncated/genome.tsv"
                }
            ],
            "output_providers": [
                "supergsl.core.output.ASTPrintOutputProvider"
            ],
            "plugins": [
                "supergsl.core.parts.provider",
                "supergsl.plugins.chopchop.plugin"
            ]
        }
