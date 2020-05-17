class Node(object):
    def child_nodes(self):
        return []


class SymbolRepository(object):

    def __init__(self):
        self._symbols = {}

    def register(self, name, table_obj):
        if name in self._symbols:
            raise Exception('Symbol table collision. Table %s is already present in global symbol table.' % name)

        self._symbols[name] = table_obj


class Program(Node):
    def __init__(self, imports, assembly_list):
        self.assembly_list = assembly_list
        self.imports = imports

    def eval(self):
        return {
            'node': 'Program',
            'imports': [
                impor.eval()
                for impor in self.imports
            ],
            'assemblies': [
                assembly.eval()
                for assembly in self.assembly_list
            ]
        }

    def child_nodes(self):
        return self.imports + self.assembly_list


class AssemblyList(Node):
    def __init__(self, assemblies):
        self.assemblies = assemblies

    def eval(self):
        return [
            assembly.eval()
            for assembly in self.assemblies
        ]

    def child_nodes(self):
        return self.assemblies


class Assembly(Node):
    def __init__(self, parts):
        self.parts = parts

    def eval(self):
        return {
            'node': 'Assembly',
            'parts': [
                part.eval()
                for part in self.parts
            ]
        }

    def child_nodes(self):
        return self.parts


class Part(Node):
    def __init__(self, identifier, slice=None):
        self.identifier = identifier
        self.operator_prefix = identifier[0]
        self.part_name = identifier[1:]
        self.slice = slice

        self.validate_part_prefix()

    def validate_part_prefix(self):
        """Validate the part prefix.

        From the GSL Paper, valid part prefixes are the following:
        g prefix gene locus gADH1
        p prefix promoter part pERG10
        t prefix terminator part tERG10
        u prefix upstream part uHO
        d prefix downstream part dHO
        o prefix open reading frame oERG10
        f prefix fusible ORF, no stop codon fERG10
        m prefix mRNA (ORF + terminator)
        """
        if self.operator_prefix not in 'gptudofm':
            raise Exception('Invalid part prefix "%s" in "%s".' % (
                self.operator_prefix,
                self.identifier
            ))

    def eval(self):
        result = {
            'node': 'Part',
            'operator_prefix': self.operator_prefix,
            'part_name': self.part_name,
            'part': self.part
        }

        if self.slice:
            result['slice'] = self.slice.eval()

        return result

    def child_nodes(self):
        if self.slice:
            return [self.slice]
        return []


class Slice(Node):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def eval(self):
        return {
            'node': 'Slice',
            'start': self.start,
            'end': self.end
        }


class ProgramImport(Node):
    def __init__(self, module_path, import_identifiers):
        self.module = module_path
        self.imports = import_identifiers

    def eval(self):
        return {
            'node': 'Import',
            'module': self.module,
            'imports': [
                progam_import.eval()
                for progam_import in self.imports
            ]
        }

    def child_nodes(self):
        return self.imports


class ProgramImportIdentifier(Node):
    def __init__(self, identifier):
        self.identifier = identifier

    def eval(self):
        return {
            'identifier': self.identifier
        }
