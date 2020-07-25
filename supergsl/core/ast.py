from __future__ import annotations
from typing import Dict, List, Optional

class Node(object):
    def child_nodes(self) -> List[Node]:
        return []


class SymbolRepository(object):

    def __init__(self):
        self._symbols = {}

    def register(self, name : str, table_obj):
        if name in self._symbols:
            raise Exception('Symbol table collision. Table %s is already present in global symbol table.' % name)

        self._symbols[name] = table_obj


class Part(Node):
    def __init__(self, identifier : str, slice : Optional[None] = None):
        self.identifier = identifier
        self.operator_prefix = identifier[0]
        self.part_name = identifier[1:]
        self.slice = slice

        self.part_type = self.get_part_type()

    def get_part_type(self) -> str:
        """Validate the part prefix.

        https://github.com/Amyris/GslCore/blob/b738b3e107b91ed50a573b48d0dcf1be69c4ce6a/src/GslCore/CommonTypes.fs#L60

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

        PART_TYPES = {
            'g': 'gene',
            'p': 'promoter',
            't': 'terminator',
            'u': 'upstream',
            'd': 'downstream',
            'o': 'orf',
            'f': 'fusible_orf',
            'm': 'mRNA'
        }

        try:
            return PART_TYPES[self.operator_prefix]
        except KeyError:
            raise Exception('Invalid part prefix "%s" in "%s".' % (
                self.operator_prefix,
                self.identifier
            ))

    def eval(self) -> dict:
        result = {
            'node': 'Part',
            'operator_prefix': self.operator_prefix,
            'part_name': self.part_name,
            #'source_part': self.source_part,
            'part_type': self.part_type
        }

        if self.slice:
            result['slice'] = self.slice.eval()

        return result

    def child_nodes(self):
        if self.slice:
            return [self.slice]
        return []


class Slice(Node):
    def __init__(self, start : int, end : int):
        self.start = start
        self.end = end

    def eval(self) -> dict:
        return {
            'node': 'Slice',
            'start': self.start,
            'end': self.end
        }


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
    def __init__(self, parts : List[Part]):
        self.parts = parts

    def eval(self) -> Dict:
        return {
            'node': 'Assembly',
            'parts': [
                part.eval()
                for part in self.parts
            ]
        }

    def child_nodes(self):
        return self.parts


class Program(Node):
    def __init__(self, imports, assembly_list : List[Assembly]):
        self.assembly_list = assembly_list
        self.imports = imports

    def eval(self) -> dict:
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

    def child_nodes(self) -> List[Node]:
        return self.imports + self.assembly_list
