from typing import TextIO
from supergsl.core.function import SuperGSLFunction
from supergsl.core.assembly import (
    AssemblyResultOutputFunction,
    AssemblyResultSet
)

from sbol2 import (
    ComponentDefinition,
    Config,
    Document,
    Sequence,
    setHomespace,

)


class SBOLOutput(AssemblyResultOutputFunction):
    """Generate SBOL document containing the assemblies."""

    name = 'sbol'


    def sanitize_identifier(self, identifier):
        """Sanitize SuperGSL Identifiers to conform to identifiers in the SBOL spec."""
        bad_chars = '[]~:'
        for c in bad_chars:
            identifier = identifier.replace(c, '_')
        return identifier

    def output(self, assemblies : AssemblyResultSet, file_handle : TextIO):
        """Output Assembly results as an SBOL document."""

        setHomespace('http://sbols.org/SuperGSL_Example/')
        Config.setOption('sbol_compliant_uris', True)
        Config.setOption('sbol_typed_uris', True)

        self.sbol_doc = Document()
        self.assembly_count = 0

        for assembly in assemblies:
            self.handle_assembly(assembly)

        file_handle.write(self.sbol_doc.writeString())


    def handle_assembly(self, assembly):
        """Add each assembly to the SBOL Document."""
        self.assembly_count += 1

        label = assembly.identifier
        if not label:
            label = 'Assembly%05d' % self.assembly_count

        assembly_comp_def = ComponentDefinition(label)
        self.sbol_doc.addComponentDefinition(assembly_comp_def)

        part_components = []
        for part in assembly.parts:
            sanitized_ident = self.sanitize_identifier(part.identifier)

            part_component = ComponentDefinition(sanitized_ident)
            part_component.roles = part.roles
            part_component.sequence = Sequence(sanitized_ident, str(part.sequence))

            part_components.append(part_component)

        assembly_comp_def.assemblePrimaryStructure(part_components)
        assembly_comp_def.compile()
