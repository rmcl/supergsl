from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunctionDeclaration

from .fuse import FusionAssembler
from .print import SuperGSLTypeDetailFunction, SuperGSLTypeHelpFunction
from .output.json_output import JSONOutput
from .output.sbol_output import SBOLOutput
from .output.genbank_output import GenBankOutput

class BuiltinAssemblersPlugin(SuperGSLPlugin):
    """Plugin stub to help register basic Assemblers."""

    def register(self, compiler_settings):
        """Register built in assemblers."""
        self.register_function('builtin', 'fuse', SuperGSLFunctionDeclaration(
            FusionAssembler, compiler_settings))

        self.register_function('builtin', 'detail', SuperGSLFunctionDeclaration(
            SuperGSLTypeDetailFunction, compiler_settings))

        self.register_function('builtin', 'help', SuperGSLFunctionDeclaration(
            SuperGSLTypeHelpFunction, compiler_settings))

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
