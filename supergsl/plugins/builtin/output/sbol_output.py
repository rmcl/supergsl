from typing import TextIO
from supergsl.core.types.part import Part
from supergsl.core.types.assembly import Assembly
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

    def output(self, assemblies : AssemblyResultSet, file_handle : TextIO) -> None:
        """Output Assembly results as an SBOL document."""

        setHomespace('http://sbols.org/SuperGSL/')
        Config.setOption('sbol_compliant_uris', True)
        Config.setOption('sbol_typed_uris', True)

        sbol_doc = Document()

        assembly_count = 1
        for assembly in assemblies:
            self.handle_assembly(sbol_doc, assembly, assembly_count)
            assembly_count += 1

        file_handle.write(sbol_doc.writeString())


    def handle_assembly(
        self,
        sbol_doc : Document,
        assembly : Assembly,
        assembly_count : int
    ) -> None:
        """Add an assembly to a SBOL Document."""

        label = assembly.identifier
        if not label:
            label = 'Assembly%05d' % assembly_count

        assembly_comp_def = ComponentDefinition(label)
        sbol_doc.addComponentDefinition(assembly_comp_def)

        part_components = []
        for part in assembly.reagents:

            if not isinstance(part, Part):
                # Exclude all reagents that are not parts
                continue

            sanitized_ident = self.sanitize_identifier(part.identifier)

            part_component = ComponentDefinition(sanitized_ident)
            part_component.roles = part.roles
            part_component.sequence = Sequence(sanitized_ident, str(part.sequence))

            part_components.append(part_component)

        assembly_comp_def.assemblePrimaryStructure(part_components)
        assembly_comp_def.compile()

    def sanitize_identifier(self, identifier : str) -> str:
        """Sanitize SuperGSL Identifiers to conform to identifiers in the SBOL spec."""
        bad_chars = '[]~:'
        for c in bad_chars:
            identifier = identifier.replace(c, '_')
        return identifier
