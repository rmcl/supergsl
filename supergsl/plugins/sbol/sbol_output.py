from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunction, SuperGSLFunctionDeclaration
from supergsl.core.assembly import AssemblyList

from sbol2 import (
    ComponentDefinition,
    Config,
    Document,
    Sequence,
    setHomespace,

)


class SBOLOutput(SuperGSLFunction):
    """Generate SBOL document containing the assemblies."""

    name = 'sbol'

    def get_arguments(self):
        return [
            ('filename', str),
            ('assemblies', list)
        ]

    def get_return_type(self):
        return type(None)

    def sanitize_identifier(self, identifier):
        """Sanitize SuperGSL Identifiers to conform to identifiers in the SBOL spec."""
        bad_chars = '[]~:'
        for c in bad_chars:
            identifier = identifier.replace(c, '_')
        return identifier

    def handle_assembly(self, assembly):
        """Add each assembly to the SBOL Document."""
        self.assembly_count += 1

        label = node.label
        if not label:
            label = 'Assembly%05d' % self.assembly_count

        assembly = ComponentDefinition(label)
        self.sbol_doc.addComponentDefinition(assembly)

        part_components = []
        for part_node in node.parts:
            part = part_node.part

            sanitized_ident = self.sanitize_identifier(part.identifier)

            part_component = ComponentDefinition(sanitized_ident)
            part_component.roles = part.roles
            part_component.sequence = Sequence(sanitized_ident, str(part.get_sequence().seq))

            part_components.append(part_component)

        assembly.assemblePrimaryStructure(part_components)
        assembly.compile()

    def execute(self, params : dict):
        """Initialize the SBOL Document."""

        setHomespace('http://sbols.org/SuperGSL_Example/')
        Config.setOption('sbol_compliant_uris', True)
        Config.setOption('sbol_typed_uris', True)

        self.sbol_doc = Document()
        self.assembly_count = 0

        assembly_list : AssemblyList = params[0]
        for assembly in assembly_list:
            self.handle_assembly(assembly)

        self.sbol_doc.write('output_sbol.xml')


class SBOLPlugin(SuperGSLPlugin):
    """Plugin stub to register SBOL related functions."""

    def register(self, compiler_settings):
        """Register SBOL functions."""
        self.register_function(
            'sbol',
            'save_sbol',
            SuperGSLFunctionDeclaration(SBOLOutput, compiler_settings))
