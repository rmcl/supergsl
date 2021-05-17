"""Entrypoint for the `sgsl-util` command used to access utility functions of the SuperGSL compiler."""
import argparse
from supergsl.core.config import load_settings
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.parts.provider import PartProviderPlugin

def handle_part_command(args):
    settings = load_settings()
    symbol_table = SymbolTable()
    part_provider = PartProviderPlugin()
    part_provider.register(symbol_table, settings)

    if args.action == 'list':
        part_provider = symbol_table.get_plugin_provider(args.part_provider)
        parts = part_provider.list_parts()

        print('\t'.join([
            'Identifier', 'Alternative Name', 'Description'
        ]))
        for part in parts:
            print('\t'.join([
                part.name,
                ','.join(part.alternative_names),
                part.description
            ]))

    else:
        raise Exception('unknown command %s' % args.action)

def main():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands',
        help='additional help',
        dest='subcommand',
        required=True)

    parser_parts = subparsers.add_parser('part', help='Interact with superGSL part providers.')
    parser_parts.add_argument('part_provider', type=str, help='The provider you wish to introspect.')
    parser_parts.add_argument('action', choices=['list'], help='The action you would like to perform')

    parser_assembly = subparsers.add_parser('assembly', help='Interact with superGSL assembly providers.')
    parser_assembly.add_argument('assembly_provider', type=str, help='The assembly provider you wish to introspect.')

    args = parser.parse_args()
    print(args)

    subcommands = {
        'part': handle_part_command
    }

    handler = subcommands[args.subcommand]
    handler(args)


if __name__ == "__main__":
    main()
