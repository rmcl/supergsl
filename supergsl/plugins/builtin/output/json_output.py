"""Implement a SuperGSL function to output Assemblies as JSON."""
from typing import Dict, List, Tuple, Type, Any, TextIO
import json
from supergsl.core.assembly import AssemblyResultOutputFunction
from supergsl.core.types.assembly import AssemblyResultSet

class JSONOutput(AssemblyResultOutputFunction):
    """Generate JSON document containing a list of parts and assemblies.

    Expected contents of JSON document:

        parts
        assemblies

    """

    def output(self, assemblies : AssemblyResultSet, file_handle : TextIO):
        parts = set()
        json_output : Dict[str, Any] = {
            'parts': [],
            'assemblies': []
        }

        for assembly_idx, assembly in enumerate(assemblies):
            assembly_parts = [
                part.identifier
                for part in assembly.parts
            ]

            json_output['assemblies'].append({
                'identifier': assembly.identifier,
                'parts': assembly_parts,
                'sequence': str(assembly.sequence)
            })

            for part in assembly.parts:
                parts.add(part)

        for part in parts:
            part_details = {
                'identifier': part.identifier,
                'sequence': str(part.sequence),
            }

            if part.has_primers:
                primers = part.get_extraction_primers()
                part_details['forward_primer'] = str(primers.forward.sequence)
                part_details['reverse_primer'] = str(primers.reverse.sequence)

            json_output['parts'].append(part_details)

        file_handle.write(json.dumps(json_output, sort_keys=True, indent=4))
