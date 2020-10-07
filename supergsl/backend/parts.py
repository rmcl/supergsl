from collections import namedtuple
from supergsl.utils import import_class
from supergsl.core.exception import ConfigurationException
from supergsl.core.config import settings
from supergsl.core.backend import BreadthFirstNodeFilteredPass


class PartSymbolTable(object):

    def __init__(self):
        self._parts = {}

        self._initialize_providers()

    def get_part(self, part_alias):
        try:
            return self._parts[part_alias]
        except KeyError:
            raise Exception('Part "%s" has not been defined.' % part_alias)

    def resolve_part(self, provider_name, part_name, alias=None):
        print('Resolving Part: %s, %s, %s' % (provider_name, part_name, alias))
        if not alias:
            alias = part_name

        if alias in self._parts:
            raise Exception('Part "%s" is already defined.' % alias)

        provider = self.resolve_provider(provider_name)
        part = provider.get_part(part_name)

        self._parts[alias] = part

    def _initialize_providers(self):
        self._providers = {}

        if 'part_providers' not in settings:
            raise ConfigurationException('No part providers have been defined. Check your supergGSL settings.')

        for provider_config in settings['part_providers']:
            print('Initializing "%s"' % provider_config['name'])
            provider_class = import_class(provider_config['provider_class'])
            provider_inst = provider_class(provider_config['name'], provider_config)

            provider_name = provider_inst.get_provider_name()
            if not provider_name:
                raise ConfigurationException('Provider "%s" does not specify a name.' % provider_class)

            self._providers[provider_name] = provider_inst

    def resolve_provider(self, provider_name):
        try:
            return self._providers[provider_name]
        except KeyError:
            raise Exception('Unknown part provider "%s".' % provider_name)


class Part(object):
    def __init__(self, name, sequence, parent_part=None):
        self.name = name
        self.sequence = sequence
        self.parent_part = parent_part

    def get_child_part(self, part_slice, new_part_name):
        """Retrieve a subsequence of this part. For example the promoter region."""
        # Obtain sliced sequence.
        left_offset = part_slice.left.x
        if part_slice.left.rel_to == 'ThreePrime':
            left_offset *= -1

        right_offset = part_slice.right.x
        if part_slice.left.rel_to == 'ThreePrime':
            left_offset *= -1

        sub_sequence = self.sequence[left_offset:right_offset]

        return Part(new_part_name, sub_sequence, self)

    def __eq__(self, other):
        return self.sequence == other.sequence

    def __hash__(self):
        return hash(self.sequence)

    def __repr__(self):
        return self.name


class PartProvider(object):
    name = None

    def __init__(self, name):
        self.name = name

    def get_provider_name(self):
        return self.name

    def get_part(self, identifier):
        """Retrieve a part from the provider.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """
        raise NotImplemented('Subclass to implement.')


class ResolvePartPass(BreadthFirstNodeFilteredPass):

    def get_node_handlers(self):
        return {
            'ProgramImport': self.visit_import_node,
            'Part': self.visit_part_node,
        }

    def before_pass(self, ast):
        self.part_symbol_table = PartSymbolTable()
        ast.symbol_registry.register('parts', self.part_symbol_table)
        return ast

    def visit_import_node(self, node):
        for program_import in node.imports:
            self.part_symbol_table.resolve_part(
                '.'.join(node.module),
                program_import.identifier
            )

    def visit_part_node(self, node):
        node.source_part = self.part_symbol_table.get_part(node.part_name)


SequencePosition = namedtuple('SequencePosition', ['x', 'rel_to', 'approximate'])

def get_promoter_len() -> int:
    """Get the configured length of a promoter region.

    ### TODO: REFACTOR THIS METHOD TO COME FROM CONFIG
    """

    return 500

def get_terminator_len() -> int:
    return 500

def get_flank_len() -> int:
    """Get the configured length of a flanking region.

    ### TODO: REFACTOR THIS METHOD TO COME FROM CONFIG
    """
    return 500


class DNASlice(object):
    def __init__(self, source_part, part_type, left, right):
        self.source_part = source_part
        self.part_type = part_type
        self.left = left
        self.right = right

    def __str__(self):
        return '%s %s %s %s' % (
            str(self.source_part),
            self.part_type,
            self.left,
            self.right
        )


class SliceAndBuildPartSequencePass(BreadthFirstNodeFilteredPass):
    """Visit each assembly part, slice and determine final sequence.

    Parts can be sliced to return subsequences of supplied genetic parts.

    Classic *fGSL* Slice Syntax:

    `<prefix>PARTNAME[<slice>]`
    1. Use the <prefix> to determine the part type slice, for example `pGAL1`.
    2. Use the <slice> to determine the sub-region of the part type for example pGAL1[1:300] will
        return the first 300 bp of the promoter region.


    Hierarchical Part Syntax

    ```
    PARTNAME[<subcomponent>]
    PARTNAME[<subcomponent>]
    ```

    Use the bracket notation to access child parts. This syntax supports infinite child components, though in practice
    more than one or two is likely to rather confusing.

    ```
    PARTNAME[<subcomponent>][<child-of-subcomponent]
    ```

    You cannot utilize *fGSL* part prefix with hierarchical parts, but you can prepend `.promoter`, `.terminator`, etc
    to access these regions based on standard GSL semantics. i.e `promoter` standards for first 500 bp by default or is
    overriden by PROMOTER_LENGTH setting.

    ```
    PARTNAME[<subcomponent>].promoter
    ```

    """

    def get_node_handlers(self):
        return {
            'Part': self.visit_part_node,
        }

    def visit_part_node(self, node):
        part = node.source_part
        part_type = node.get_part_type()

        child_part = self.slice_part_by_part_type(part, part_type)

        """
        # derive a new part from the previous part based on slice notation
        if node.part_slice:
            print('SLICE', node.part_slice)
            raise NotImplemented('Slice not implemented yet.')
        """


    def slice_part_by_part_type(self, part, part_type):
        part_type_slice = self.build_part_type_slice(
            part,
            part_type
        )

        new_part_name = part_type + part.name
        new_part = part.get_child_part(part_type_slice, new_part_name)


    def build_part_type_slice(self, part, part_type):
        """Build the slice of a part based on the requested part type.

        parts often have a part type specified by the prefix, for example
        p for promoter.

        Refer to translateGenePrefix in GslCore for reference logic
        https://github.com/Amyris/GslCore/blob/d2c613907d33b110a2f53021146342234e0d8f3b/src/GslCore/DnaCreation.fs#L53

        """
        if part_type == 'promoter':
            return DNASlice(
                source_part=part,
                part_type='promoter',
                left=SequencePosition(
                    x=get_promoter_len(),
                    rel_to='FivePrime',
                    approximate=True
                ),
                right=SequencePosition(
                    x=-1,
                    rel_to='FivePrime',
                    approximate=False
                )
            )

        elif part_type == 'upstream':
            return DNASlice(
                source_part=part,
                part_type='promoter',
                left=SequencePosition(
                    x=(-1 * get_flank_len()),
                    rel_to='FivePrime',
                    approximate=True
                ),
                right=SequencePosition(
                    x=-1,
                    rel_to='FivePrime',
                    approximate=False
                )
            )

        elif part_type == 'downstream':
            return DNASlice(
                source_part=part,
                part_type='promoter',
                left=SequencePosition(
                    x=0,
                    rel_to='ThreePrime',
                    approximate=False
                ),
                right=SequencePosition(
                    x=get_flank_len(),
                    rel_to='ThreePrime',
                    approximate=True
                )
            )

        elif part_type == 'gene':
            return DNASlice(
                source_part=part,
                part_type='promoter',
                left=SequencePosition(
                    x=0,
                    rel_to='FivePrime',
                    approximate=False
                ),
                right=SequencePosition(
                    x=-1,
                    rel_to='ThreePrime',
                    approximate=False
                )
            )

        elif part_type == 'terminator':
            return DNASlice(
                source_part=part,
                part_type='promoter',
                left=SequencePosition(
                    x=0,
                    rel_to='ThreePrime',
                    approximate=False
                ),
                right=SequencePosition(
                    x=get_terminator_len(),
                    rel_to='ThreePrime',
                    approximate=True
                )
            )

        else:
            raise Exception('"%s" not implemented yet.' % part_type)
