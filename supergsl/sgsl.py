"""Entrypoint for the `sgsl` command used to invoke the superGSL compiler."""
import argparse
from supergsl.core.config import load_settings
from supergsl.core.pipeline import CompilerPipeline
from supergsl.core.output import OutputPipeline
from supergsl.repl import SuperGSLShell



def main():
    settings = load_settings()
    compiler_pipeline = CompilerPipeline(settings)
    output_pipeline = OutputPipeline(settings)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file",
        help="The input source code file to process",
        type=str,
        default=None,
        nargs='?')

    parser.add_argument(
        "-f",
        "--output-format",
        action='append',
        help="Specify compiler output formats",
        type=str,
        required=False)

    args = parser.parse_args()

    print(args.input_file)
    if not args.input_file:
        SuperGSLShell().start()
    else:
        output_pipeline.validate_args(args)

        print('Compiling "%s".' % args.input_file)

        with open(args.input_file, 'r') as input_file_fp:
            result = compiler_pipeline.compile(input_file_fp.read())

        print('Compiling Complete. Storing Output.')

        output_pipeline.run(result, args)


if __name__ == "__main__":
    main()
