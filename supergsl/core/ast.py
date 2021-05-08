from __future__ import annotations
from typing import cast, Dict, List, Optional, Any, Union, Type, Tuple

from supergsl.core.types import SuperGSLType, Collection
from supergsl.core.function import SuperGSLFunction, SuperGSLFunctionDeclaration

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


class SymbolReference(Node):
    def __init__(self, identifier : str, slice : Optional[Slice], invert : bool):
        self.identifier = identifier
        self.slice = slice
        self.invert = invert

        self.referenced_object : Optional[SuperGSLType] = None

    def set_table_reference(self, symbol_table, identifier):
        """Set the SuperGSL object that has been associated with this symbol."""
        self.symbol_table = symbol_table
        self.symbol_table_identifier = identifier

    def eval(self) -> SuperGSLType:
        symbol = self.symbol_table.lookup(self.symbol_table_identifier)
        return symbol.eval(self)

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
    def __init__(self, identifier: str, type_declaration : Optional[TypeDeclaration], value : ListDeclaration):
        self.identifier : str = identifier
        self.type_declaration : Optional[TypeDeclaration] = type_declaration
        self.value : ListDeclaration = value

    def eval(self) -> Tuple[str, SuperGSLType]:
        result = self.value.eval()

        if self.type_declaration:
            # Todo: Implement type declaration check
            pass

        return self.identifier, result

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
    def __init__(self, items : List[SuperGSLType]):
        self.items = items

    def eval(self) -> Collection:
        return Collection([
            part.eval()
            for part in self.items
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
    def __init__(self, parts : List[SymbolReference], label : Optional[str] = None):
        self.parts : List[SymbolReference] = parts
        self.label : Optional[str] = label

    def eval(self) -> List[SymbolReference]:
        return [
            part.eval()
            for part in self.parts
        ]

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
    def __init__(
        self, identifier : str,
        children : DefinitionList,
        params : Optional[List[str]] = None,
        label : str = None
    ):
        self.identifier = identifier
        self.children = children
        self.params = params
        self.label = label

        self.function_declaration : Optional[SuperGSLFunctionDeclaration] = None

    def set_function_declaration(self, function_def : SuperGSLFunctionDeclaration):
        self.function_declaration = function_def

    def eval(self) -> SuperGSLFunction:
        if not self.function_declaration:
            raise Exception('Function has not been defined.')

        function_inst = self.function_declaration.eval(self)

        expected_return_type = function_inst.get_return_type()
        print('expected', expected_return_type)

        eval_params = {
            'children': []
        }

        if self.params:
            ## TODO: THIS IS NOT RIGHT YET!
            print('params', self.params)
            for idx in range(len(self.params)):
                print(self.params[idx])
                eval_params[idx] = self.params[idx].eval()

        if self.children:
            eval_params['children'] = [
                child.eval()
                for child in self.children.definitions
            ]

        function_result = function_inst.execute(eval_params)

        if not isinstance(function_result, expected_return_type):
            raise Exception('Unexpected function return type. Expected: %s, Actual: %s' % (
                type(function_result),
                expected_return_type))

        return function_result

    def to_dict(self) -> Dict:
        results = {
            'node': 'FunctionInvocation',
            'identifier': self.identifier,
            'children': self.children.to_dict() if self.children else None,
            'label': self.label
        }

        if self.params:
            results['params'] = [
                param.to_dict()
                for param in self.params
            ]
        return results

    def get_definition_list(self):
        return self.children

    def child_nodes(self):
        all_child_nodes = []
        if self.children:
            all_child_nodes.append(self.children)
        if self.params:
            all_child_nodes.extend(self.params)
        return all_child_nodes

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
