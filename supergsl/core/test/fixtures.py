ex1 = (
    '''uHO ; pADH1''',
    'IDENTIFIER',
    'SEMICOLON',
    'IDENTIFIER'
)
ex2 = (
    '''uHO ; pADH1 ; gERG10[1:728] ; dHO'''
    'IDENTIFIER',
    'SEMICOLON',
    'IDENTIFIER',
    'LEFT_BRACKET'
    'NUMBER',
    'COLON',
    'NUMBER',
    'RIGHT_BRACKET',
    'SEMICOLON',
)

# I guess ### means marker based on a pragma. We're not doing this for now.
#ex4 = """uHO ; pADH1 ; gERG10[1:728] ; ### dHO"""
ex3 = """
from S288C import ADHA, ERG10, HO
uHO ; pADH1 ; gERG10[1:728] ; dHO
"""

ex4 = """
from S288C import ADHA, ERG10, HO
uHO ; pADH1 ; gERG10[1:728] ; ### ; dHO
"""
