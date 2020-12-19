from __future__ import annotations
from typing import cast, Dict, List, Optional, Any, Union


class Node(object):
    def child_nodes(self) -> List[Node]:
        return []

    def eval(self) -> Dict[str, Any]:
        return {}

    def replace_child_node(self, cur_node, new_node):
        raise NotImplementedError('"%s" does not implement `replace_child_node`.' % self)


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


class SlicePosition(Node):
    def __init__(self, index: int, postfix : str, approximate : boolean):
        self.index = index
        self.postfix = postfix
        self.approximate = approximate

    def eval(self) -> Dict[str, Any]:
        return {
            'node': 'SlicePosition',
            'index': self.index,
            'postfix': self.postfix,
            'approximate': self.approximate
        }

class Slice(Node):
    def __init__(self, start : int, end : int):
        self.start = start
        self.end = end

    def eval(self) -> Dict[str, Any]:
        return {
            'node': 'Slice',
            'start': self.start.eval(),
            'end': self.end.eval()
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

    def child_nodes(self) -> List[Node]:
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

    def child_nodes(self) -> List[Node]:
        return cast(List[Node], self.imports)


Definition = Union[
    'Assembly',
    'FunctionInvocation'
]

class DefinitionList(Node):
    def __init__(self, definitions: List[Definition]):
        self.definitions = definitions

    def add_definition(self, definition : Definition):
        self.definitions.append(definition)

    def eval(self) -> dict:
        return {
            'node': 'DefinitionList',
            'items': [
                definition.eval()
                for definition in self.definitions
            ]
        }

    def child_nodes(self) -> List[Node]:
        return cast(List[Node], self.definitions)

    def replace_child_node(self, old_node, new_node):
        """Replace an existing child node with a new replacement node."""
        self.definitions = [
            new_node if node is old_node else node
            for node in self.definitions
        ]

class Assembly(Node):
    def __init__(self, parts : List[Part], label : Optional[str] = None):
        self.parts : List[Part] = parts
        self.label : Optional[str] = label

    def eval(self) -> Dict:
        return {
            'node': 'Assembly',
            'parts': [
                part.eval()
                for part in self.parts
            ],
            'label': self.label
        }

    def child_nodes(self) -> List[Node]:
        return cast(List[Node], self.parts)


class FunctionInvocation(Node):
    def __init__(self, identifier : str, children : DefinitionList, params : Optional[List[str]] = None, label : str = None):
        self.identifier = identifier
        self.children = children
        self.params = params
        self.label = label

    def eval(self) -> Dict:
        return {
            'node': 'FunctionInvocation',
            'identifier': self.identifier,
            'children': self.children.eval() if self.children else None,
            'params': self.params,
            'label': self.label
        }

    def get_definition_list(self):
        return self.child_definition_list

    def child_nodes(self):
        if self.child_definition_list:
            return [self.child_definition_list]
        return []


class NucleotideConstant(Node):
    def __init__(self, sequence : str):
        self.sequence = sequence

    def eval(self) -> dict:
        return {
            'node': 'NucleotideConstant',
            'sequence': self.sequence
        }


class Program(Node):
    def __init__(self, imports : List[ProgramImport], definitions : DefinitionList):
        self.definitions : DefinitionList = definitions
        self.imports : List[ProgramImport] = imports

    def eval(self) -> dict:
        return {
            'node': 'Program',
            'imports': [
                impor.eval()
                for impor in self.imports
            ],
            'definitions': self.definitions.eval()
        }

    def child_nodes(self) -> List[Node]:
        children : List[Node] = cast(List[Node], self.imports.copy())
        children.append(cast(Node, self.definitions))
        return children
