"""Implement a SuperGSL function to output Assemblies as JSON."""
from typing import Dict, List, Tuple, Type, Any
import json
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunction, SuperGSLFunctionDeclaration
from supergsl.core.types.assembly import AssemblyResultSet

class JSONOutput(SuperGSLFunction):
    """Generate JSON document containing a list of parts and assemblies.

    Expected contents of JSON document:

        parts
        assemblies

    """

    return_type = None

    @classmethod
    def get_arguments(cls) -> List[Tuple[str, Type]]:
        return [
            ('assemblies', AssemblyResultSet)
            #('filename', str),
        ]

    def execute(self, params : dict):
        parts = set()
        json_output : Dict[str, Any] = {
            'parts': [],
            'assemblies': []
        }

        assembly_list : AssemblyResultSet = params['assemblies']
        for assembly_idx, assembly in enumerate(assembly_list):
            assembly_sequence = assembly.get_sequence()
            assembly_parts = [
                part.identifier
                for part in assembly.get_required_parts()
            ]

            json_output['assemblies'].append({
                'identifier': assembly_idx,
                'parts': assembly_parts,
                'sequence': str(assembly_sequence)
            })

            for part in assembly.get_required_parts():
                parts.add(part)

        for part in parts:
            part_details = {
                'name': part.identifier,
                'sequence': str(part.get_sequence().seq),
            }

            if part.has_primers:
                primers = part.get_extraction_primers()
                primers.forward.get_sequence()
                part_details['forward_primer'] = primers.forward.get_sequence()
                part_details['reverse_primer'] = primers.reverse.get_sequence()

            json_output['parts'].append(part_details)

        print(json.dumps(json_output, sort_keys=True, indent=4))



class JSONOutputPlugin(SuperGSLPlugin):
    """Plugin stub to help register basic Assemblers."""

    def register(self, compiler_settings : dict):
        """Register built in assemblers."""
        self.register_function(
            'json',
            'output_json',
            SuperGSLFunctionDeclaration(JSONOutput, compiler_settings))
