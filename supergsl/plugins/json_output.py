from supergsl.core.backend import BreadthFirstNodeFilteredPass


class JSONOutputPass(BreadthFirstNodeFilteredPass):
    """Generate JSON document containing a list of parts and assemblies.

    Expected contents of JSON document:
    
        parts
        assemblies

    """

    def get_node_handlers(self):
        return {
            'Assembly': self.visit_assembly_node,
        }

    def before_pass(self, ast):
        """Initialize the SBOL Document."""
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

    def visit_assembly_node(self, node):

        assembly_parts = [
            part.identifier
            for part in node.parts
        ]

        assembly_sequence = ''.join([
            str(part.source_part.sequence.seq)
            for part in node.parts
        ])

        assembly_idx = len(self.json_output['assemblies'])
        self.json_output['assemblies'].append({
            'identifier': assembly_idx,
            'parts': assembly_parts,
            'sequence': assembly_sequence
        })
