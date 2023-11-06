from typing import List, Optional


class SuperGSLIntegrationFixtures(object):
    """Fixtures to assist in execution of integration tests."""

    def get_supergsl_settings(self, extra_plugins: Optional[List[str]] = None) -> dict:
        """Return a basic nested dictionary of settings.

        extra_plugins (list): a list of paths to plugins to include in the settings.
        """

        settings : dict = {
            "part_providers": [
                {
                    "name": "truncated.S288C",
                    "provider_class": "supergsl.plugins.file_part.FeatureTableWithFastaPartProvider",
                    "fasta_file_path": "supergsl/tests/fixtures/s288c_truncated/genome.fa.gz",
                    "feature_file_path": "supergsl/tests/fixtures/s288c_truncated/genome.tsv"
                }
            ],
            "output_providers": [
                "supergsl.core.output.ASTPrintOutputProvider"
            ]
        }
        settings['plugins'] = [
            "supergsl.core.parts.provider",
            "supergsl.plugins.chopchop.plugin"
        ] + (extra_plugins or [])

        return settings
