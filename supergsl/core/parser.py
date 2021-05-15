"""Define a parser of the SuperGSL language."""
from rply import ParserGenerator
from supergsl.core.constants import (
    UNAMBIGUOUS_DNA_SEQUENCE,
    UNAMBIGUOUS_PROTEIN_SEQUENCE
)

from . import ast
from .exception import ParsingError




class ParserState(object):
    def __init__(self, filename=None):
        self.filename = filename

        # Flag set to True when at least one assembly has been
        # defined. If we reach the end of compiler input and this
        # is false then we should raise an error.
        self.one_assembly_defined = False


class ParserBuilder(object):
    # A list of all token names accepted by the parser.
    ACCEPTED_TOKENS = (
        'FROM',
        'IMPORT',
        'AS',
        'LET',

        'OPEN_CURLY_BRACKET',
        'CLOSE_CURLY_BRACKET',

        'FORWARD_SLASH',

        'OPEN_PAREN',
        'CLOSE_PAREN',

        'OPEN_BRACKET',
        'CLOSE_BRACKET',
        'COLON',
        'SEMICOLON',
        'COMMA',
        'PERIOD',
        'NUMBER',
        'IDENTIFIER',
        'EQUAL',
        'TILDE',
        'EXCLAMATION',
        'DOLLAR_SIGN',
    )

    def __init__(self):
        self.pg = ParserGenerator(self.ACCEPTED_TOKENS)
        self.build_parser()

    def build_parser(self):
        """Define the parser rules."""

        # rply has it's own style which does not conform to pylint's expectations.
        # pylint: disable=W0613,C0103,W0612,C0301

        @self.pg.production('program : import_list definition_list')
        @self.pg.production('program : definition_list')
        def program(state, p):
            imports = []
            if len(p) == 2:
                # we have at least one import
                imports = p[0]
                definition_list = p[1]
            else:
                definition_list = p[0]

            return ast.Program(imports, definition_list)

        @self.pg.production('import_list : import_list import')
        @self.pg.production('import_list : import')
        def program_import_list(state, p):
            if len(p) == 2:
                p[0].append(p[1])
                return p[0]

            return [p[0]]

        @self.pg.production('import : FROM import_module IMPORT import_identifiers')
        def program_import(state, p):
            return ast.Import(p[1], p[3])

        @self.pg.production('import_module : import_module PERIOD IDENTIFIER')
        @self.pg.production('import_module : IDENTIFIER')
        def program_import_module(state, p):
            if len(p) == 3:
                return p[0] + [p[2].value]

            return [p[0].value]

        @self.pg.production('import_identifiers : import_identifiers COMMA import_identifier')
        @self.pg.production('import_identifiers : import_identifier')
        def program_import_identifiers(state, p):
            if len(p) == 3:
                p[0].append(p[2])
                return p[0]

            return [p[0]]

        @self.pg.production('import_identifier : IDENTIFIER AS IDENTIFIER')
        @self.pg.production('import_identifier : IDENTIFIER')
        def import_identifier(state, p):
            import_identifier = p[0].value
            if len(p) == 3:
                alias_identifier = p[2].value
            else:
                alias_identifier = None

            return ast.ImportIdentifier(import_identifier, alias_identifier)

        @self.pg.production('definition_list : definition_list definition')
        @self.pg.production('definition_list : definition')
        def definition_list(state, p):
            if len(p) == 2:
                p[0].add_definition(p[1])
                return p[0]

            return ast.DefinitionList([p[0]])

        @self.pg.production('definition : variable_definition')
        @self.pg.production('definition : function_invoke')
        @self.pg.production('definition : assembly')
        def definition(state, p):
            return p[0]

        @self.pg.production('variable_definition : LET IDENTIFIER EQUAL symbol_reference')
        @self.pg.production('variable_definition : LET IDENTIFIER EQUAL function_invoke')
        @self.pg.production('variable_definition : LET IDENTIFIER EQUAL list_declaration')
        @self.pg.production('variable_definition : LET IDENTIFIER type_declaration EQUAL list_declaration')
        def variable_definition(state, p):
            variable_identifier = p[1].value
            if len(p) == 4:
                variable_type = None
                definition = p[3]
            elif len(p) == 5:
                variable_type = p[2]
                definition = p[4]

            return ast.VariableDeclaration(variable_identifier, variable_type, definition)

        @self.pg.production('type_declaration : OPEN_BRACKET IDENTIFIER CLOSE_BRACKET')
        def type_declaration(state, p):
            return ast.TypeDeclaration(p[1].value)

        @self.pg.production('list_declaration : OPEN_BRACKET list_items CLOSE_BRACKET')
        def list_declaration(state, p):
            return ast.ListDeclaration(p[1])

        @self.pg.production('list_items : list_item')
        @self.pg.production('list_items : list_item COMMA list_items')
        def list_items(state, p):
            new_list = [p[0]]
            if len(p) == 3:
                new_list.extend(p[2])
            return new_list

        @self.pg.production('list_item : symbol_reference')
        @self.pg.production('list_item : nucleotide_constant')
        @self.pg.production('list_item : amino_acid_constant')
        def list_item(state, p):
            return p[0]

        @self.pg.production('function_name_and_label : IDENTIFIER')
        @self.pg.production('function_name_and_label : IDENTIFIER COLON IDENTIFIER')
        def function_name_and_label(state, p):
            if len(p) == 3:
                # return Name and label
                return p[2].value, p[0].value

            return p[0].value, None

        @self.pg.production('function_invoke : function_name_and_label function_parameter_block')
        @self.pg.production('function_invoke : function_name_and_label OPEN_CURLY_BRACKET definition_list CLOSE_CURLY_BRACKET')
        def function_invoke(state, p):
            if len(p) == 2:
                return ast.FunctionInvocation(p[0][0], None, p[1], p[0][1])

            return ast.FunctionInvocation(p[0][0], p[2], None, p[0][1])


        @self.pg.production('function_parameter_block : OPEN_PAREN function_parameters CLOSE_PAREN')
        @self.pg.production('function_parameter_block : OPEN_PAREN CLOSE_PAREN')
        def function_param_block(state, p):
            if len(p) == 2:
                return None

            return p[1]

        @self.pg.production('function_parameters : symbol_reference')
        @self.pg.production('function_parameters : symbol_reference COMMA function_parameters')
        def function_parameters(state, p):
            x = [p[0]]
            if len(p) == 3:
                x.extend(p[2])
            return x

        @self.pg.production('assembly : assembly_list')
        @self.pg.production('assembly : IDENTIFIER COLON assembly_list')
        def assembly(state, p):
            if len(p) == 1:
                return ast.Assembly(p[0])
            elif len(p) == 3:
                return ast.Assembly(p[2], label=p[0].value)

        @self.pg.production('assembly_list : assembly_list SEMICOLON assembly_list_item')
        @self.pg.production('assembly_list : assembly_list_item')
        def assembly_list(state, p):
            if len(p) == 1:
                return [p[0]]
            elif len(p) == 3:
                p[0].append(p[2])
                return p[0]

        @self.pg.production('assembly_list_item : symbol_reference')
        @self.pg.production('assembly_list_item : nucleotide_constant')
        def assembly_list_item(state, p):
            return p[0]

        @self.pg.production('nucleotide_constant : FORWARD_SLASH IDENTIFIER FORWARD_SLASH')
        def nucleotide_constant(state, p):
            return ast.SequenceConstant(p[1].value, UNAMBIGUOUS_DNA_SEQUENCE)

        @self.pg.production('amino_acid_constant : FORWARD_SLASH DOLLAR_SIGN IDENTIFIER FORWARD_SLASH')
        def protein_constant(state, p):
            return ast.SequenceConstant(p[1].value, UNAMBIGUOUS_PROTEIN_SEQUENCE)

        @self.pg.production('symbol_reference : symbol_identifier OPEN_BRACKET slice_index CLOSE_BRACKET')
        @self.pg.production('symbol_reference : symbol_identifier')
        def symbol_reference(state, p):
            identifier, invert = p[0]
            part_slice = None
            if len(p) == 4:
                part_slice = p[2]

            return ast.SymbolReference(identifier, part_slice, invert)

        @self.pg.production('symbol_identifier : EXCLAMATION IDENTIFIER')
        @self.pg.production('symbol_identifier : IDENTIFIER')
        def symbol_identifier(state, p):
            if len(p) == 2:
                invert = True
                identifier = p[1].value
            else:
                invert = False
                identifier = p[0].value

            return (identifier, invert)

        @self.pg.production('slice_index : slice_position COLON slice_position')
        def index_slice(state, p):
            return ast.Slice(p[0], p[2])

        @self.pg.production('slice_position : TILDE slice_coordinates')
        @self.pg.production('slice_position : slice_coordinates')
        def slice_position(state, p):
            if len(p) == 2:
                approximate = True
                position_index, postfix = p[1]
            else:
                approximate = False
                position_index, postfix = p[0]

            return ast.SlicePosition(position_index, postfix, approximate)

        @self.pg.production('slice_coordinates : NUMBER')
        @self.pg.production('slice_coordinates : NUMBER IDENTIFIER')
        def slice_coordinates(state, p):
            position_index = int(p[0].value)
            postfix = None
            if len(p) == 2:
                postfix = p[1].value
                if postfix not in ['S', 'E']:
                    raise ParsingError('Slice postfix can only be "S" or "E"')

            return (position_index, postfix)

        @self.pg.error
        def error_handle(state, lookahead):
            reached_end_of_file = lookahead.value == '$end'
            if reached_end_of_file:
                if not state.one_assembly_defined:
                    raise ParsingError('At least one assembly must be defined.')

            raise ParsingError(
                'An error occurred parsing source document at %s' % lookahead.source_pos)

    def parse(self, tokens):
        parser = self.pg.build()

        parser_state = ParserState()
        return parser.parse(tokens, state=parser_state)
