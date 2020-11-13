from rply import ParserGenerator
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

        'OPEN_CURLY_BRACKET',
        'CLOSE_CURLY_BRACKET',

        'OPEN_PAREN',
        'CLOSE_PAREN',

        'OPEN_BRACKET',
        'CLOSE_BRACKET',
        'COLON',
        'SEMICOLON',
        'COMMA',
        'PERIOD',
        'NUMBER',
        'IDENTIFIER'
    )

    def __init__(self):
        self.pg = ParserGenerator(self.ACCEPTED_TOKENS)
        self.build_parser()

    def build_parser(self):
        """Define the parser rules."""
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
            else:
                return [p[0]]

        @self.pg.production('import : FROM import_module IMPORT import_identifiers')
        def program_import(state, p):
            return ast.ProgramImport(p[1], p[3])

        @self.pg.production('import_module : import_module PERIOD IDENTIFIER')
        @self.pg.production('import_module : IDENTIFIER')
        def program_import_module(state, p):
            if len(p) == 3:
                return p[0] + [p[2].value]
            else:
                return [p[0].value]

        @self.pg.production('import_identifiers : import_identifiers COMMA IDENTIFIER')
        @self.pg.production('import_identifiers : IDENTIFIER')
        def program_import_identifiers(state, p):
            if len(p) == 3:
                pi = ast.ProgramImportIdentifier(p[2].value)
                p[0].append(pi)
                return p[0]
            else:
                #print(p[0].source_pos)
                return [ast.ProgramImportIdentifier(p[0].value)]


        @self.pg.production('definition_list : definition_list function_invoke')
        @self.pg.production('definition_list : definition_list assembly')
        @self.pg.production('definition_list : function_invoke')
        @self.pg.production('definition_list : assembly')
        def definition_list(state, p):
            if len(p) == 2:
                p[0].extend(p[1])
                return p[0]
            else:
                return [p[0]]


        @self.pg.production('function_invoke : IDENTIFIER OPEN_PAREN CLOSE_PAREN OPEN_CURLY_BRACKET definition_list CLOSE_CURLY_BRACKET')
        @self.pg.production('function_invoke : IDENTIFIER OPEN_CURLY_BRACKET definition_list CLOSE_CURLY_BRACKET')
        def function_invoke(state, p):
            # Functions can be called with or without parameters enclosed by parenthesis
            if len(p) == 4:
                return ast.FunctionInvocation(p[0].value, p[2])
            elif len(p) == 6:
                raise Exception('AHHHHHHH')


        @self.pg.production('assembly : part_list')
        @self.pg.production('assembly : IDENTIFIER COLON part_list')
        def assembly(state, p):
            if len(p) == 1:
                return ast.Assembly(p[0])
            elif len(p) == 3:
                return ast.Assembly(p[2], label=p[0].value)

        @self.pg.production('part_list : part_list SEMICOLON part')
        @self.pg.production('part_list : part')
        def part_list(state, p):
            if len(p) == 1:
                return [p[0]]
            elif len(p) == 3:
                p[0].append(p[2])
                return p[0]

            assert('Cant reach this point.')

        @self.pg.production('part : IDENTIFIER slice')
        def sliced_part(state, p):
            return ast.Part(p[0].value, p[1])

        @self.pg.production('part : IDENTIFIER')
        def simple_part(state, p):
            return ast.Part(p[0].value)

        @self.pg.production('slice : OPEN_BRACKET index_slice CLOSE_BRACKET')
        def slice(state, p):
            return p[1]

        @self.pg.production('index_slice : NUMBER COLON NUMBER')
        def index_slice(state, p):
            return ast.Slice(p[0].value, p[2].value)

        @self.pg.error
        def error_handle(state, lookahead):
            reached_end_of_file = lookahead.value == '$end'
            if reached_end_of_file:
                if not state.one_assembly_defined:
                    raise ParsingError('At least one assembly must be defined.')

            raise ValueError(state)

    def parse(self, tokens):
        parser = self.pg.build()

        parser_state = ParserState()
        return parser.parse(tokens, state=parser_state)
