from lexer import Lexer
from parser import Parser

ex1 = """
uHO ; pADH1
"""
ex3 = """
uHO ; pADH1 ; gERG10[1:728] ; dHO"""

# I guess ### means marker based on a pragma. We're not doing this for now.
#ex4 = """uHO ; pADH1 ; gERG10[1:728] ; ### dHO"""
ex6 = """
from S288C import ADHA, ERG10, HO
uHO ; pADH1 ; gERG10[1:728] ; dHO
"""

ex2 = """
from S288C import ADHA, ERG10, HO
uHO ; pADH1 ; gERG10[1:728] ; ### dHO
"""

lexer = Lexer().get_lexer()
tokens = lexer.lex(ex6)

pg = Parser()
pg.parse()

parser = pg.get_parser()
result = parser.parse(tokens).eval()

import pprint
pprint.pprint(result)
