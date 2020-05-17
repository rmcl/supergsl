class Node(object):
    def child_nodes(self):
        return []


class Root(object):

    def __init__(self):
        self._symbols = {}

    def add_symbol_table(self, name, table_obj):
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
        self.slice = slice

    def eval(self):
        result = {
            'node': 'Part',
            'identifier': self.identifier
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
        self.import_identifiers = import_identifiers

    def eval(self):
        return {
            'node': 'Import',
            'module': self.module,
            'imports': [
                import_identifier.eval()
                for import_identifier in self.import_identifiers
            ]
        }

    def child_nodes(self):
        return self.import_identifiers


class ProgramImportIdentifier(Node):
    def __init__(self, identifier):
        self.identifier = identifier

    def eval(self):
        return {
            'identifier': self.identifier
        }
