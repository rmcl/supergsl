from supergsl.core.output import OutputProvider


class JSONOutputPass(OutputProvider):
    """Generate JSON document containing a list of parts and assemblies.

    Expected contents of JSON document:

        parts
        assemblies

    """

    name = 'json'

    def get_node_handlers(self):
        return {
            'Assembly': self.visit_assembly_node,
            'Part': self.visit_part_node,
        }

    def before_pass(self, ast):
        """Initialize the SBOL Document."""
        self.parts = set()
        self.json_output = {
            'parts': [],
            'assemblies': []
        }
        return ast

    def after_pass(self, ast):
        import pprint
        pprint.pprint(self.json_output)
        #self.sbol_doc.write('output_sbol.xml')
        return ast

    def visit_part_node(self, node):
        part = node.part
        if part in self.parts:
            return

        self.parts.add(part)
        self.json_output['parts'].append({
            'name': part.name,
            'slice_of_parent': str(part.slice_of_parent[0]),
            'forward_primer': str(part.forward_primer.seq),
            'reverse_primer': str(part.reverse_primer.seq),
            'sequence': str(part.sequence.seq),
        })

    def visit_assembly_node(self, node):
        assembly_parts = [
            part.identifier
            for part in node.parts
        ]

        assembly_sequence = ''.join([
            str(part.sequence.seq)
            for part in node.parts
        ])

        assembly_idx = len(self.json_output['assemblies'])
        self.json_output['assemblies'].append({
            'identifier': assembly_idx,
            'parts': assembly_parts,
            'sequence': assembly_sequence
        })
