import sys
import textwrap
from io import StringIO
from unittest import TestCase
from supergsl.utils import get_logo, display_symbol_table
from supergsl.core.symbol_table import SymbolTable

class UtilitiesTestCases(TestCase):
    """Test the helper methods."""

    def test_get_logo(self):
        logo_data = get_logo()
        self.assertEquals(len(logo_data.split('\n')), 9)

    def test_display_symbol_table(self):
        st = SymbolTable('GLOBAL', None)
        st.insert('BOB', 'HELLO')
        st.nested_scope('BOOM')
        st.insert('BOB', 'TAMBORINE')

        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            display_symbol_table(st)
            result = sys.stdout.getvalue()
        except Exception:
            sys.stdout = old_stdout
            raise
        finally:
            sys.stdout = old_stdout

        self.assertEqual(result, textwrap.dedent("""
            Table: GLOBAL ----
                Key: BOB, Type: <class 'str'>

                Table: BOOM ----
                --- End Table: BOOM ----
            --- End Table: GLOBAL ----
            """))
