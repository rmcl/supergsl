from __future__ import annotations
from typing import (
    Dict,
    List,
    Optional,
    Any,
    Union,
    Type,
    cast,
)

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


class SymbolReference(Node):
    def __init__(self, identifier : str, slice : Optional[Slice], invert : bool):
        self.identifier = identifier
        self.slice = slice
        self.invert = invert

    def to_dict(self) -> Dict[str, Any]:
        return {
            'node': 'SymbolReference',
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
    def __init__(self, module_path : List[str], import_identifiers : List[ImportIdentifier]):
        self.module_path = module_path
        self.imports = import_identifiers

    def to_dict(self):
        return {
            'node': 'Import',
            'module': self.module_path,
            'imports': [
                progam_import.to_dict()
                for progam_import in self.imports
            ]
        }

    def child_nodes(self) -> List[Node]:
        return cast(List[Node], self.imports)

    def get_node_label(self):
        return '%s:%s' % (self.__class__.__name__, '.'.join(self.module_path))


Definition = Union[
    'Assembly',
    'VariableDeclaration',
    # Todo: Maybe re-enable FunctionInvoke as a definition here
    # 'FunctionInvocation'
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
    def __init__(
        self,
        identifier: str,
        type_declaration : Optional[TypeDeclaration],
        value : ListDeclaration
    ):
        self.identifier : str = identifier
        self.type_declaration : Optional[TypeDeclaration] = type_declaration
        self.value : ListDeclaration = value

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
    def __init__(self, items : List[Node]):
        self.item_nodes : List[Node] = items

    def to_dict(self) -> dict:
        return {
            'node': 'ListDeclaration',
            'items': [
                item_node.to_dict()
                for item_node in self.item_nodes
            ]
        }

    def child_nodes(self) -> List[Node]:
        return self.item_nodes


class Assembly(Node):
    def __init__(self, parts : List[SymbolReference], label : Optional[str] = None):
        self.symbol_references : List[SymbolReference] = parts
        self.label : Optional[str] = label

    def to_dict(self) -> Dict:
        return {
            'node': 'Assembly',
            'parts': [
                symbol_reference.to_dict()
                for symbol_reference in self.symbol_references
            ],
            'label': self.label
        }

    def child_nodes(self) -> List[Node]:
        return cast(List[Node], self.symbol_references)


class FunctionInvocation(Node):
    """AST node representing function calls."""
    def __init__(
        self,
        identifier : str,
        children : DefinitionList,
        params : List[Any],
        label : Optional[str]
    ):
        self.identifier = identifier
        self.children = children
        self.params = params
        self.label = label

    def to_dict(self) -> Dict:
        results : Dict[str, Any] = {
            'node': 'FunctionInvocation',
            'identifier': self.identifier,
            'children': self.children.to_dict() if self.children else None,
            'params': None,
            'label': self.label
        }

        if self.params:
            results['params'] = [
                param.to_dict()
                for param in self.params
            ]
        return results

    def child_nodes(self):
        all_child_nodes = []
        if self.children:
            all_child_nodes.append(self.children)
        if self.params:
            all_child_nodes.extend(self.params)
        return all_child_nodes

class Constant(Node):
    def __init__(self, value : str, constant_type : Type):
        self.value = value
        self.constant_type = constant_type

    def to_dict(self) -> dict:
        return {
            'node': 'Constant',
            'type': self.constant_type,
            'value': self.value
        }

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
    def __init__(self, imports : List[Import], definitions : Optional[DefinitionList]):
        self.definitions : Optional[DefinitionList] = definitions
        self.imports : List[Import] = imports

    def to_dict(self) -> dict:
        return {
            'node': 'Program',
            'imports': [
                impor.to_dict()
                for impor in self.imports
            ],
            'definitions': self.definitions.to_dict() if self.definitions else None
        }

    def child_nodes(self) -> List[Node]:
        children : List[Node] = cast(List[Node], self.imports.copy())
        if self.definitions:
            children.append(cast(Node, self.definitions))
        return children
