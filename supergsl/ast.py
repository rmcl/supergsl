class Program(object):
    def __init__(self, imports, assembly_list):
        self.assembly_list = assembly_list
        self.imports = imports

    def eval(self):
        return {
            'node': 'Program',
            'imports': self.imports.eval() if self.imports else None,
            'assemblies': self.assembly_list.eval()
        }


class AssemblyList(object):
    def __init__(self, assemblies):
        self.assemblies = assemblies

    def eval(self):
        return [
            assembly.eval()
            for assembly in self.assemblies
        ]


class Assembly(object):
    def __init__(self, part_list):
        self.part_list = part_list

    def eval(self):
        return {
            'node': 'Assembly',
            'parts': self.part_list.eval()
        }

class PartList(object):

    def __init__(self, parts):
        self.parts = parts

    def append(self, part):
        self.parts.append(part)

    def eval(self):
        return [
            part.eval()
            for part in self.parts
        ]


class Part(object):
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

class Slice(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def eval(self):
        return {
            'start': self.start,
            'end': self.end
        }

class ProgramImportList(object):
    def __init__(self, program_imports):
        self.program_imports = program_imports

    def append(self, program_import):
        self.program_imports.append(program_import)

    def eval(self):
        return [
            program_import.eval()
            for program_import in self.program_imports
        ]

class ProgramImport(object):
    def __init__(self, module_path, import_identifiers):
        self.module = module_path
        self.import_identifiers = import_identifiers

    def eval(self):
        return {
            'node': 'Import',
            'module': self.module,
            'imports': self.import_identifiers.eval()
        }

class ProgramImportIdentifiers(object):
    def __init__(self, identifiers):
        self.identifiers = identifiers

    def append(self, identifier):
        self.identifiers.append(identifier)

    def eval(self):
        return [
            identifier
            for identifier in self.identifiers
        ]
