from supergsl.core.output import OutputProvider
from sbol import (
    setHomespace,
    Config,
    Document,
    ComponentDefinition,
    Sequence,

    SO_GENE,
    SO_PROMOTER,
    SO_TERMINATOR,
    SO_MISC,
    SO_CDS,
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

        self.part_type_role_map = {
            'gene': [SO_GENE],
            'promoter': [SO_PROMOTER],
            'terminator': [SO_TERMINATOR],
            'upstream': [SO_MISC],
            'downstream': [SO_MISC],
            'orf': ['0000236'], # SO_ORF seems not to be defined in libsbol. but SO is 0000236 http://www.sequenceontology.org/browser/current_release/term/SO:0000236
            'fusible_orf': [SO_CDS],
            'mRNA': [],
        }

        setHomespace('http://sbols.org/SUPERGSL_Example/')
        Config.setOption('sbol_compliant_uris', True)
        Config.setOption('sbol_typed_uris', True)

        self.sbol_doc = Document()

        return ast

    def after_pass(self, ast):
        self.sbol_doc.write('output_sbol.xml')
        return ast

    def visit_assembly_node(self, node):

        assembly = ComponentDefinition('Assembly')
        self.sbol_doc.addComponentDefinition(assembly)

        part_components = []
        for part in node.parts:
            part_component = ComponentDefinition(part.identifier)
            part_component.roles = self.part_type_role_map[part.get_part_type()]
            part_component.sequence = Sequence(part.identifier, str(part.source_part.sequence.seq))

            part_components.append(part_component)

        assembly.assemblePrimaryStructure(part_components)
        assembly.compile()
