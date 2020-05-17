from rply import ParserGenerator
from . import ast


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
                'NUMBER',
                'IDENTIFIER'
            ]
        )

    def parse(self):

        @self.pg.production('program : import_list assembly')
        @self.pg.production('program : assembly')
        def program(p):
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
        def program_import_list(p):
            if len(p) == 2:
                p[0].append(p[1])
                return p[0]
            else:
                return [p[0]]

        @self.pg.production('import : FROM IDENTIFIER IMPORT import_identifiers')
        def program_import(p):
            return ast.ProgramImport(p[1].value, p[3])

        @self.pg.production('import_identifiers : import_identifiers COMMA IDENTIFIER')
        @self.pg.production('import_identifiers : IDENTIFIER')
        def program_import_identifiers(p):
            if len(p) == 3:
                pi = ast.ProgramImportIdentifier(p[2].value)
                p[0].append(pi)
                return p[0]
            else:
                return [ast.ProgramImportIdentifier(p[0].value)]

        @self.pg.production('assembly : part_list')
        def assembly(p):
            return ast.Assembly(p[0])

        @self.pg.production('part_list : part_list SEMICOLON part')
        @self.pg.production('part_list : part')
        def part_list(p):
            if len(p) == 1:
                return [p[0]]
            elif len(p) == 3:
                p[0].append(p[2])
                return p[0]

            assert('Cant reach this point.')

        @self.pg.production('part : IDENTIFIER slice')
        def sliced_part(p):
            return ast.Part(p[0].value, p[1])

        @self.pg.production('part : IDENTIFIER')
        def simple_part(p):
            return ast.Part(p[0].value)

        @self.pg.production('slice : OPEN_BRACKET NUMBER COLON NUMBER CLOSE_BRACKET')
        def slice(p):
            return ast.Slice(p[1].value, p[3].value)

        @self.pg.error
        def error_handle(token):
            raise ValueError(token)

    def get_parser(self):
        return self.pg.build()
