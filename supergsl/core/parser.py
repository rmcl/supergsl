from rply import ParserGenerator
from . import ast


class ParserState(object):
    def __init__(self):
        pass


class ParserBuilder(object):
    def __init__(self):
        self.pg = ParserGenerator(
            # A list of all token names accepted by the parser.
            [
                'FROM',
                'IMPORT',

                'OPEN_BRACKET',
                'CLOSE_BRACKET',
                'COLON',
                'SEMICOLON',
                'COMMA',
                'PERIOD',
                'NUMBER',
                'IDENTIFIER'
            ]
        )

        self.build_parser()

    def build_parser(self):

        @self.pg.production('program : import_list assembly')
        @self.pg.production('program : assembly')
        def program(state, p):
            imports = None
            if len(p) == 2:
                # we have at least one import
                imports = p[0]
                assembly = [p[1]]
            else:
                assembly = [p[0]]

            return ast.Program(imports, assembly)

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
                return [ast.ProgramImportIdentifier(p[0].value)]

        @self.pg.production('assembly : part_list')
        def assembly(state, p):
            return ast.Assembly(p[0])

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
            return ast.Slice(p[1].value, p[3].value)

        @self.pg.error
        def error_handle(token):
            raise ValueError(token)

    def parse(self, tokens):
        parser = self.pg.build()

        parser_state = ParserState()

        print(tokens)
        return parser.parse(tokens, state=parser_state)
