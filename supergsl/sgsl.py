"""Entrypoint for the `sgsl` command used to invoke the superGSL compiler."""
import argparse
from supergsl.core.config import load_settings
from supergsl.core.pipeline import CompilerPipeline
from supergsl.core.exception import SuperGSLError
from supergsl.repl import SuperGSLShell
import pprint


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file",
        help="The input source code file to process",
        type=str,
        default=None,
        nargs='?')

    parser.add_argument(
        "-s", "--start-shell-on-error",
        help="If an error occurs during execution of SuperGSL program then start the repl shell.",
        default=False,
        action='store_true')

    args = parser.parse_args()

    compiler_settings = load_settings()
    compiler_pipeline = CompilerPipeline(compiler_settings)

    if not args.input_file:
        SuperGSLShell(compiler_pipeline).start()
    else:
        print('Compiling "%s".' % args.input_file)

        with open(args.input_file, 'r') as input_file_fp:
            source_code = input_file_fp.read()

        try:
            compiler_pipeline.compile(source_code)
        except SuperGSLError as error:
            if args.start_shell_on_error:
                SuperGSLShell(compiler_pipeline).start()
            else:
                raise error


    print('Compiling Complete.')

if __name__ == "__main__":
    main()
