from __future__ import annotations
from typing import cast, Dict, List, Optional, Any

class Node(object):
    def child_nodes(self) -> List[Node]:
        return []

    def eval(self) -> Dict[str, Any]:
        return {}


class SymbolRepository(object):

    def __init__(self):
        self._symbols = {}

    def register(self, name : str, table_obj):
        if name in self._symbols:
            raise Exception('Symbol table collision. Table %s is already present in global symbol table.' % name)

        self._symbols[name] = table_obj

    def get_table(self, name : str):
        try:
            return self._symbols[name]
        except KeyError:
            raise Exception('Unknown symbol table "%s".' % name)


class Slice(Node):
    def __init__(self, start : int, end : int):
        self.start = start
        self.end = end

    def eval(self) -> Dict[str, Any]:
        return {
            'node': 'Slice',
            'start': self.start,
            'end': self.end
        }


class Part(Node):
    def __init__(self, identifier : str, slice : Optional[Slice] = None):
        self.identifier = identifier
        self.slice = slice

    def eval(self) -> Dict[str, Any]:
        result : Dict[str, Any] = {
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

    def __str__(self):
        return '(%s) %s' % (self.identifier, self.part)


class ProgramImportIdentifier(Node):
    def __init__(self, identifier : str):
        self.identifier : str = identifier

    def eval(self) -> dict:
        return {
            'identifier': self.identifier
        }


class ProgramImport(Node):
    def __init__(self, module_path : str, import_identifiers : List[ProgramImportIdentifier]):
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


class Assembly(Node):
    def __init__(self, parts : List[Part], label : Optional[str] = None):
        self.parts = parts
        self.label = label

    def eval(self) -> Dict:
        return {
            'node': 'Assembly',
            'parts': [
                part.eval()
                for part in self.parts
            ],
            'label': self.label
        }

    def child_nodes(self):
        return self.parts


class AssemblyBlock(Node):
    def __init__(self, assembly_type, assemblies: List[Assembly]):
        self.assembly_type = assembly_type
        self.assemblies = assemblies
        print(self.assemblies)

    def eval(self) -> Dict:
        return {
            'node': 'AssemblyBlock',
            'assembly_type': self.assembly_type,
            'assemblies': [
                assembly.eval()
                for assembly in self.assemblies
            ]
        }

    def child_nodes(self):
        return self.assemblies

class Program(Node):
    def __init__(self, imports : List[ProgramImport], assembly_blocks : List[AssemblyBlock]):
        self.assembly_blocks = assembly_blocks
        self.imports = imports

    def eval(self) -> dict:
        return {
            'node': 'Program',
            'imports': [
                impor.eval()
                for impor in self.imports
            ],
            'assembly_blocks': [
                assembly_block.eval()
                for assembly_block in self.assembly_blocks
            ]
        }

    def child_nodes(self) -> List[Node]:
        return cast(List[Node], self.imports) + cast(List[Node], self.assembly_blocks)
