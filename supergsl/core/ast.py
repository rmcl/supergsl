from __future__ import annotations
from typing import cast, Dict, List, Optional, Any, Union

from supergsl.core.parts.part import Part as CorePart
from supergsl.core.types import SuperGSLType, Collection

# rply has it's own style which does not conform to pylint's expectations.
# pylint: disable=E1136

class Node(object):
    def child_nodes(self) -> List['Node']:
        return []

    def to_dict(self) -> Dict[str, Any]:
        return {}

    def replace_child_node(self, cur_node, new_node):
        raise NotImplementedError('"%s" does not implement `replace_child_node`.' % self)

    def get_node_label(self):
        return str(self.__class__.__name__)

class SlicePosition(Node):
    def __init__(self, index: int, postfix : str, approximate : bool):
        self.index = index
        self.postfix = postfix
        self.approximate = approximate

    def to_dict(self) -> Dict[str, Any]:
        return {
            'node': 'SlicePosition',
            'index': self.index,
            'postfix': self.postfix,
            'approximate': self.approximate
        }

    def get_slice_pos_str(self):
        return '%s%d%s' % (
            '~' if self.approximate else '',
            self.index,
            self.postfix if self.postfix else ''
        )

class Slice(Node):
    def __init__(self, start : SlicePosition, end : SlicePosition):
        self.start = start
        self.end = end

    def to_dict(self) -> Dict[str, Any]:
        return {
            'node': 'Slice',
            'start': self.start.to_dict(),
            'end': self.end.to_dict()
        }

    def get_slice_str(self):
        return '%s:%s' % (
            self.start.get_slice_pos_str(),
            self.end.get_slice_pos_str()
        )


class Part(Node):
    def __init__(self, identifier : str, slice : Optional[Slice], invert : bool):
        self.identifier = identifier
        self.slice = slice
        self.invert = invert

        self.part : Optional[CorePart] = None
        self.parent_parts : List[CorePart] = []

    def eval(self) -> SuperGSLType:
        return self.part

    def to_dict(self) -> Dict[str, Any]:
        return {
            'node': 'Part',
            'identifier': self.identifier,
            'invert': self.invert,
            'slice': self.slice.to_dict() if self.slice else None
        }

    def child_nodes(self) -> List[Node]:
        if self.slice:
            return [self.slice]
        return []

    def __str__(self):
        return self.identifier

    def get_node_label(self):
        return '%s:%s' % (self.__class__.__name__, self.identifier)


class ImportIdentifier(Node):
    def __init__(self, identifier : str, alias : Optional[str]):
        self.identifier : str = identifier
        self.alias : Optional[str] = alias

    def to_dict(self) -> dict:
        return {
            'identifier': self.identifier,
            'alias': self.alias
        }

    def get_node_label(self):
        label = '%s: %s' % (
            self.__class__.__name__,
            self.identifier,
        )
        if self.alias:
            label = label + ' (%s)' % self.alias
        return label

class Import(Node):
    def __init__(self, module_path : str, import_identifiers : List[ImportIdentifier]):
        self.module = module_path
        self.imports = import_identifiers

    def to_dict(self):
        return {
            'node': 'Import',
            'module': self.module,
            'imports': [
                progam_import.to_dict()
                for progam_import in self.imports
            ]
        }

    def child_nodes(self) -> List[Node]:
        return cast(List[Node], self.imports)

    def get_node_label(self):
        return '%s:%s' % (self.__class__.__name__, '.'.join(self.module))


Definition = Union[
    'Assembly',
    'FunctionInvocation'
]

class DefinitionList(Node):
    def __init__(self, definitions: List[Definition]):
        self.definitions = definitions

    def add_definition(self, definition : Definition):
        self.definitions.append(definition)

    def to_dict(self) -> dict:
        return {
            'node': 'DefinitionList',
            'items': [
                definition.to_dict()
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

class TypeDeclaration(Node):
    def __init__(self, identifier: str):
        self.identifier : str = identifier

    def to_dict(self) -> dict:
        return {
            'node': 'TypeDeclaration',
            'identifier': self.identifier
        }

class VariableDeclaration(Node):
    def __init__(self, identifier: str, type_declaration : Optional[TypeDeclaration], value : ListDeclaration):
        self.identifier : str = identifier
        self.type_declaration : Optional[TypeDeclaration] = type_declaration
        self.value : ListDeclaration = value

    def eval(self):
        return AssignmentOp(self.identifier, self.value.eval())

    def to_dict(self) -> dict:
        return {
            'node': 'VariableDeclaration',
            'identifier': self.identifier,
            'value': self.value.to_dict(),
            'type_declaration': (
                self.type_declaration.to_dict() if self.type_declaration else None
            )
        }

    def child_nodes(self) -> List[Node]:
        return [cast(Node, self.value)]


class ListDeclaration(Node):
    def __init__(self, items : List[Part]):
        self.items = items

    def eval(self) -> SuperGSLType:
        if len(self.items) == 0:
            raise Exception('Empty list not supported.')

        return Collection([
            item.eval()
            for item in self.items
        ])

    def to_dict(self) -> dict:
        return {
            'node': 'ListDeclaration',
            'items': [
                part.to_dict()
                for part in self.items
            ]
        }

    def child_nodes(self) -> List[Node]:
        return cast(List[Node], self.items)

class Assembly(Node):
    def __init__(self, parts : List[Part], label : Optional[str] = None):
        self.parts : List[Part] = parts
        self.label : Optional[str] = label

    def to_dict(self) -> Dict:
        return {
            'node': 'Assembly',
            'parts': [
                part.to_dict()
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

    def to_dict(self) -> Dict:
        return {
            'node': 'FunctionInvocation',
            'identifier': self.identifier,
            'children': self.children.to_dict() if self.children else None,
            'params': self.params,
            'label': self.label
        }

    def get_definition_list(self):
        return self.children

    def child_nodes(self):
        if self.children:
            return [self.children]
        else:
            return []


class SequenceConstant(Node):
    def __init__(self, sequence : str, sequence_type : str):
        self.sequence = sequence
        self.sequence_type = sequence_type

    def to_dict(self) -> dict:
        return {
            'node': 'SequenceConstant',
            'type': self.sequence_type,
            'sequence': self.sequence
        }


class Program(Node):
    def __init__(self, imports : List[Import], definitions : DefinitionList):
        self.definitions : DefinitionList = definitions
        self.imports : List[Import] = imports

    def to_dict(self) -> dict:
        return {
            'node': 'Program',
            'imports': [
                impor.to_dict()
                for impor in self.imports
            ],
            'definitions': self.definitions.to_dict()
        }

    def child_nodes(self) -> List[Node]:
        children : List[Node] = cast(List[Node], self.imports.copy())
        children.append(cast(Node, self.definitions))
        return children
