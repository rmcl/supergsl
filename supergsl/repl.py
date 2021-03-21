"""Entrypoint for the `sgsl` command used to invoke the superGSL compiler."""
from types import SimpleNamespace
import textwrap
from prompt_toolkit import PromptSession
from prompt_toolkit.history import FileHistory
from prompt_toolkit.auto_suggest import AutoSuggestFromHistory

from supergsl.core.config import load_settings
from supergsl.core.pipeline import CompilerPipeline
from supergsl.core.output import OutputPipeline
from supergsl.core.exception import SuperGSLError
from supergsl.utils import get_logo


class SuperGSLShell:
    intro = (
        get_logo() +
        'Welcome to the SuperGSL shell.   Type help or ? to list commands.\n'
    )
    prompt = '[sGSL] '

    def start(self):
        settings = load_settings()
        self.compiler_pipeline = CompilerPipeline(settings)
        self.output_pipeline = OutputPipeline(settings)
        self.last_ast = None

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
        print(textwrap.dedent(
            """
            SuperGSL Help!

            Output Providers
                .<provider-name>

            To-Exit the shell: "Ctrl-D"

            """))

    def loop(self):

        while True:
            try:
                inp = self.prompt_for_input()
            except KeyboardInterrupt:
                continue  # Control-C pressed. Try again.
            except EOFError:
                break  # Control-D pressed.

            print('PRINT COMMAND', inp)
            if inp[0] == '?':
                self.run_help()
            elif inp[0] == '.':
                args = SimpleNamespace()
                args.output_format = [inp[1:]]

                self.output_pipeline.validate_args(args)
                self.output_pipeline.run(self.last_ast, args)

            else:
                try:
                    self.last_ast = self.run_compiler(inp)
                except SuperGSLError as error:
                    print('ERROR', error)

        print("GoodBye!")
