"""Entrypoint for the `sgsl` command used to invoke the superGSL compiler."""
import argparse
from supergsl.core.config import load_settings
from supergsl.core.pipeline import CompilerPipeline
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

    args = parser.parse_args()

    compiler_settings = load_settings()
    if not args.input_file:
        SuperGSLShell(compiler_settings).start()
    else:
        print('Compiling "%s".' % args.input_file)

        compiler_pipeline = CompilerPipeline(compiler_settings)
        with open(args.input_file, 'r') as input_file_fp:
            compiler_pipeline.compile(input_file_fp.read())

    print('Compiling Complete.')

if __name__ == "__main__":
    main()
