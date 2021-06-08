"""Define a parser of the SuperGSL language."""
from rply import ParserGenerator
from supergsl.core.constants import (
    UNAMBIGUOUS_DNA_SEQUENCE,
    UNAMBIGUOUS_PROTEIN_SEQUENCE,
    NUMBER_CONSTANT
)

from . import ast
from .exception import ParsingError




class ParserState(object):
    def __init__(self, filename=None):
        self.filename = filename


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
        'AMINO_ACID_SEQUENCE',
        'EQUAL',
        'TILDE',
        'EXCLAMATION',
    )

    def __init__(self):
        self.pg = ParserGenerator(self.ACCEPTED_TOKENS)
        self.build_parser()

    def build_parser(self):
        """Define the parser rules."""

        # rply has it's own style which does not conform to pylint's expectations.
        # pylint: disable=W0613,C0103,W0612,C0301

        @self.pg.production('program : import_list definition_list')
        @self.pg.production('program : import_list')
        @self.pg.production('program : definition_list')
        def program(state, p):
            imports = []
            definition_list = None
            if len(p) == 2:
                # we have at least one import
                imports = p[0]
                definition_list = p[1]
            else:
                if isinstance(p[0], ast.DefinitionList):
                    definition_list = p[0]
                else:
                    imports = p[0]

            return ast.Program(imports, definition_list)

        @self.pg.production('import_list : import_list import')
        @self.pg.production('import_list : import')
        def program_import_list(state, p):
            if len(p) == 2:
                p[0].append(p[1])
                return p[0]

            return [p[0]]

        @self.pg.production('import : FROM import_module IMPORT OPEN_PAREN import_identifiers CLOSE_PAREN')
        @self.pg.production('import : FROM import_module IMPORT import_identifiers')
        def program_import(state, p):
            if len(p) == 6:
                import_module = p[1]
                import_identifiers = p[4]
            if len(p) == 4:
                import_module = p[1]
                import_identifiers = p[3]

            return ast.Import(import_module, import_identifiers)

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

        @self.pg.production('variable_definition : LET IDENTIFIER EQUAL definition_item')
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

        @self.pg.production('list_items : definition_item')
        @self.pg.production('list_items : definition_item COMMA list_items')
        def list_items(state, p):
            new_list = [p[0]]
            if len(p) == 3:
                new_list.extend(p[2])
            return new_list

        @self.pg.production('definition_item : symbol_reference')
        @self.pg.production('definition_item : nucleotide_constant')
        @self.pg.production('definition_item : amino_acid_constant')
        @self.pg.production('definition_item : number_constant')
        def definition_item(state, p):
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

        @self.pg.production('function_parameters : definition_item')
        @self.pg.production('function_parameters : definition_item COMMA function_parameters')
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

        @self.pg.production('assembly_list : assembly_list SEMICOLON definition_item')
        @self.pg.production('assembly_list : definition_item')
        def assembly_list(state, p):
            if len(p) == 1:
                return [p[0]]
            elif len(p) == 3:
                p[0].append(p[2])
                return p[0]

        @self.pg.production('number_constant : NUMBER')
        def number_constant(state, p):
            return ast.Constant(p[0].value, NUMBER_CONSTANT)

        @self.pg.production('nucleotide_constant : FORWARD_SLASH IDENTIFIER FORWARD_SLASH')
        def nucleotide_constant(state, p):
            return ast.SequenceConstant(p[1].value, UNAMBIGUOUS_DNA_SEQUENCE)

        @self.pg.production('amino_acid_constant : FORWARD_SLASH AMINO_ACID_SEQUENCE FORWARD_SLASH')
        def protein_constant(state, p):
            return ast.SequenceConstant(p[1].value[1:], UNAMBIGUOUS_PROTEIN_SEQUENCE)

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
            raise ParsingError(
                'An error occurred parsing source document at %s' % lookahead.source_pos)

    def parse(self, tokens):
        parser = self.pg.build()

        parser_state = ParserState()
        return parser.parse(tokens, state=parser_state)
