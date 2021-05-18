"""Provide an interactive shell for SuperGSL."""
from types import SimpleNamespace
import textwrap
from prompt_toolkit import PromptSession
from prompt_toolkit.history import FileHistory
from prompt_toolkit.auto_suggest import AutoSuggestFromHistory

from supergsl.core.config import load_settings
from supergsl.core.pipeline import CompilerPipeline
from supergsl.core.exception import SuperGSLError
from supergsl.utils import get_logo, display_symbol_table


class SuperGSLShell:
    intro = (
        get_logo() +
        'Welcome to the SuperGSL shell.   Type help or ? to list commands.\n'
    )
    prompt = '[sGSL] '

    def __init__(self, compiler_settings : dict):
        self.settings = compiler_settings

    def start(self):
        self.compiler_pipeline = CompilerPipeline(self.settings)

        self.prompt_session = PromptSession(
            auto_suggest=AutoSuggestFromHistory(),
            history=FileHistory('sgsl_history.txt'))

        print(self.intro)
        self.loop()

    def prompt_for_input(self):
        return self.prompt_session.prompt(
            self.prompt,
            complete_while_typing=True)

    def run_compiler(self, source_code):
        return self.compiler_pipeline.compile(source_code)

    def run_help(self):
        """Display a help message."""
        print(textwrap.dedent(
            """
            SuperGSL Help!

            Commands:
            .symbols: Display a the symbol table containing all loaded objects.

            To-Exit the shell: "Ctrl-D"

            """))

    def loop(self):
        """Loop awaiting user input."""
        while True:
            try:
                inp = self.prompt_for_input()
            except KeyboardInterrupt:
                continue  # Control-C pressed. Try again.
            except EOFError:
                break  # Control-D pressed.

            if inp.strip() == '':
                continue

            if inp[0] == '?':
                self.run_help()
            elif inp == '.symbols':
                display_symbol_table(self.compiler_pipeline.get_symbol_table())
            else:
                try:
                    self.run_compiler(inp)
                except SuperGSLError as error:
                    print('ERROR', error)

        print("GoodBye!")
