"""Entrypoint for the `sgsl-util` command- access utility functions of the SuperGSL compiler."""

import argparse
from supergsl.core.config import load_settings
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.plugin import PluginProvider
from supergsl.core.parts.provider import PartProviderPlugin, PartProvider
from supergsl.utils import setup_supergsl_config_files


class CommandGroup:

    def build_parser(self, main_subparsers):
        """Build argument parser specific to parts."""
        raise NotImplementedError()

    def handle_command(self, args):
        """Handle related subcommands."""


class SetupCommands(CommandGroup):
    """Commands to initailize a SuperGSL environment including configuration files.

    Establish a `.supergsl` folder in the users home directory and initialize a
    default `supergsl-config.json` as well as a `sgsl-lib/` directory to store cached
    sequences and other data.
    """

    def build_parser(self, main_subparsers):
        """Build argument parser specific to parts."""
        setup_parser = main_subparsers.add_parser(
            'setup', help='Initialize your SuperGSL environment.')
        return setup_parser

    def handle_command(self, args):
        return setup_supergsl_config_files()


class PartProviderCommands(CommandGroup):
    """Provide CLI access to Parts and their Providers"""

    def __init__(self):
        self.symbol_table = SymbolTable()
        self.settings = load_settings()
        self.plugins = PluginProvider(self.symbol_table, self.settings)

    def handle_provider_list(self, args):
        providers = self.symbol_table.get_providers()

        for provider in providers:
            if not isinstance(provider, PartProvider):
                continue

            print(provider.name)


    def handle_part_list(self, args):
        """List parts from a specific provider. Some Providers may not support this."""
        providers = self.symbol_table.get_providers_for_path(args.part_provider)

        parts = []
        for provider in providers:
            if not isinstance(provider, PartProvider):
                continue

            parts += provider.list_parts()

        print('\t'.join([
            'Identifier', 'Alternative Name', 'Description'
        ]))
        for part in parts:
            print('\t'.join([
                part.name,
                ','.join(part.alternative_names),
                part.description
            ]))

    def handle_part_detail(self, args):
        """Display details about a specific part."""
        part = self.symbol_table.resolve_symbol(
            args.part_provider,
            args.part_identifier)

        sequence = part.get_sequence()
        print(part.identifier)
        print('Sequence Length: {}'.format(len(sequence)))
        print('Description: {}'.format(part.description))
        print('Roles:')
        for role in part.roles:
            print(' * {}'.format(role))
        print('----\n{}...\n----\n'.format(sequence[0:50].seq))


    def build_parser(self, main_subparsers):
        """Build argument parser specific to parts."""

        parser_parts = main_subparsers.add_parser(
            'parts', help='Interact with superGSL part providers.')
        #parser_parts.add_argument(
        #    'part_provider',
        #    type=str,
        #    help='The provider you wish to introspect.',
        #    )

        part_action_subparsers = parser_parts.add_subparsers(
            title='part actions',
            description='valid subcommands to interact with parts',
            help='additional help for parts',
            dest='part_subcommand',
            required=True)

        part_action_subparsers.add_parser(
            'providers', help='List all  available part providers.')
        part_action_subparsers.add_parser(
            'list', help='List parts available from this provider.')

        parser_part_detail = part_action_subparsers.add_parser(
            'detail', help='Retrieve details about a specific part.')
        parser_part_detail.add_argument(
            'part_identifier', type=str, help='The identifier of the part of interest.')

        return parser_parts

    def handle_command(self, args):
        """Handle Part related subcommands."""
        if args.part_subcommand == 'providers':
            self.handle_provider_list(args)
        elif args.part_subcommand == 'list':
            self.handle_part_list(args)
        elif args.part_subcommand == 'detail':
            self.handle_part_detail(args)

        else:
            raise Exception('unknown command %s' % args.action)

def main():
    """Main entrypoint for SuperGSL CLI Utilities."""
    parser = argparse.ArgumentParser()

    part_cli = PartProviderCommands()
    setup_cli = SetupCommands()

    main_subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands',
        help='additional help',
        dest='subcommand',
        required=True)

    part_cli.build_parser(main_subparsers)
    setup_cli.build_parser(main_subparsers)

    args = parser.parse_args()
    subcommands = {
        'parts': part_cli.handle_command,
        'setup': setup_cli.handle_command
    }

    handler = subcommands[args.subcommand]
    handler(args)


if __name__ == "__main__":
    main()
