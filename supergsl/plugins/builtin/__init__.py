from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunctionDeclaration

from .fuse import FusionAssembler
from .oligos import SyntheticOligoAssembler
from .print import SuperGSLTypePrintFunction

class BuiltinAssemblersPlugin(SuperGSLPlugin):
    """Plugin stub to help register basic Assemblers."""

    def register(self, compiler_settings):
        """Register built in assemblers."""
        self.register_function('builtin', 'fuse', SuperGSLFunctionDeclaration(
            FusionAssembler, compiler_settings))

        self.register_function('builtin', 'print', SuperGSLFunctionDeclaration(
            SuperGSLTypePrintFunction, compiler_settings))

        self.register_function('synthesis', 'synthetic_oligos', SuperGSLFunctionDeclaration(
            SyntheticOligoAssembler, compiler_settings))
