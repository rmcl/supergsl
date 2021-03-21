"""Entrypoint for the `sgsl-util` command used to access utility functions of the SuperGSL compiler."""
import os
import json
from pathlib import Path
import argparse
from supergsl.core.config import load_settings
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.parts.provider import PartProviderPlugin



def handle_setup_config_files(args):
    """Setup SuperGSL config files and directories."""
    conf_dir = Path('~/.supergsl')
    conf_dir = conf_dir.expanduser()
    lib_dir = conf_dir.joinpath('sgsl-lib/')
    conf_file = conf_dir.joinpath('config.json')
    example_config_file = Path(__file__).parent.joinpath('../supergsl-config.json.example')

    try:
        Path.mkdir(conf_dir)
    except FileExistsError as error:
        print('Cannot setup config files at "%s" because that directory already exists.' % conf_dir)
        return

    Path.mkdir(lib_dir, exist_ok=True)

    json_config = json.load(open(example_config_file))
    json_config['lib_dir'] = lib_dir.as_posix()

    with open(conf_file, 'w+') as fp:
        json.dump(json_config, fp, indent=4)
        fp.write('\n')

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

    parser_parts = subparsers.add_parser('setup', help='Setup SuperGSL config files and directories.')

    parser_parts = subparsers.add_parser('part', help='Interact with superGSL part providers.')
    parser_parts.add_argument('part_provider', type=str, help='The provider you wish to introspect.')
    parser_parts.add_argument('action', choices=['list'], help='The action you would like to perform')

    parser_assembly = subparsers.add_parser('assembly', help='Interact with superGSL assembly providers.')
    parser_assembly.add_argument('assembly_provider', type=str, help='The assembly provider you wish to introspect.')

    args = parser.parse_args()
    print(args)

    subcommands = {
        'setup': handle_setup_config_files,
        'part': handle_part_command
    }

    handler = subcommands[args.subcommand]
    handler(args)


if __name__ == "__main__":
    main()
