from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunctionDeclaration

from .fuse import FusionAssembler
from .oligos import SyntheticOligoAssembler
from .print import SuperGSLTypeDetailFunction, SuperGSLTypeHelpFunction
from .declare import AssemblyDeclarationFunction
from .output.json_output import JSONOutput
from .output.sbol_output import SBOLOutput
from .output.genbank_output import GenBankOutput


class BuiltinPlugin(SuperGSLPlugin):
    """Plugin stub to help register basic Assemblers."""

    def register(self, compiler_settings):
        """Register built in assemblers."""
        self.register_function('builtin', 'fuse', SuperGSLFunctionDeclaration(
            FusionAssembler, compiler_settings))

        self.register_function('synthesis', 'synthetic_oligos', SuperGSLFunctionDeclaration(
            SyntheticOligoAssembler, compiler_settings))

        self.register_function('builtin', 'detail', SuperGSLFunctionDeclaration(
            SuperGSLTypeDetailFunction, compiler_settings))

        self.register_function('builtin', 'help', SuperGSLFunctionDeclaration(
            SuperGSLTypeHelpFunction, compiler_settings))

        self.register_function('builtin', 'declare', SuperGSLFunctionDeclaration(
            AssemblyDeclarationFunction, compiler_settings))

        self.register_function(
            'builtin',
            'output_json',
            SuperGSLFunctionDeclaration(JSONOutput, compiler_settings))

        self.register_function(
            'builtin',
            'output_sbol',
            SuperGSLFunctionDeclaration(SBOLOutput, compiler_settings))

        self.register_function(
            'builtin',
            'output_genbank',
            SuperGSLFunctionDeclaration(GenBankOutput, compiler_settings))
