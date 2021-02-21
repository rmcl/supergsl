from supergsl.core.output import OutputProvider
from sbol2 import (
    ComponentDefinition,
    Config,
    Document,
    Sequence,
    setHomespace,

)


class SBOLOutputPass(OutputProvider):
    """Generate SBOL document containing the assemblies."""

    name = 'sbol'

    def get_node_handlers(self):
        return {
            'Assembly': self.visit_assembly_node,
        }

    def before_pass(self, ast):
        """Initialize the SBOL Document."""

        setHomespace('http://sbols.org/SuperGSL_Example/')
        Config.setOption('sbol_compliant_uris', True)
        Config.setOption('sbol_typed_uris', True)

        self.sbol_doc = Document()
        self.assembly_count = 0


    def after_pass(self, ast):
        self.sbol_doc.write('output_sbol.xml')

    def sanitize_identifier(self, identifier):
        """SBOL is really particular about part names."""
        bad_chars = '[]~:'
        for c in bad_chars:
            identifier = identifier.replace(c, '_')
        return identifier

    def visit_assembly_node(self, node):
        self.assembly_count += 1

        label = node.label
        if not label:
            label = 'Assembly%05d' % self.assembly_count

        #assembly = Component(label, SBO_DNA)
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
