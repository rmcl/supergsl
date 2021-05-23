from typing import Dict, Tuple, List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import Restriction, RestrictionBatch

from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord

from supergsl.core.exception import PartNotFoundError, SuperGSLError
from supergsl.core.assembly import AssemblerBase
from supergsl.core.constants import THREE_PRIME, SO_HOMOLOGOUS_REGION
from supergsl.core.parts.provider import PartProvider
from supergsl.core.types.assembly import AssemblyDeclaration, AssemblyList
from supergsl.core.types.position import SeqPosition
from supergsl.core.types.part import Part

class InvalidBioBrickError(SuperGSLError):
    pass

# From the docs, "BioBrick RFC[10] it must not contain the following
# restriction sites, as these are unique to the prefix and suffix"
# Taken from http://parts.igem.org/Help:Standards/Assembly/RFC10
BIOBRICK_ILLEGAL_RESTRICTION_SITES = ['EcoRI', 'XbaI', 'SpeI', 'PstI', 'NotI']

BIOBRICK_PREFIX_CODING_SEQUENCE = Seq('GAATTCGCGGCCGCTTCTAG')
BIOBRICK_PREFIX_SEQUENCE = Seq('GAATTCGCGGCCGCTTCTAGAG')
BIOBRICK_SUFFIX_SEQUENCE = Seq('TACTAGTAGCGGCCGCTGCAG')

def load_biobrick_constant_sequences():
    file_path = 'supergsl/plugins/igem/biobrick.fa'
    records = SeqIO.parse(file_path, "fasta")

    return {
        record.name: record
        for record in records
    }


def check_is_valid_biobrick(sequence):
    """A valid biobrick part does not have any disallowed restriction sites except in
    the prefix and suffix.

    Return True if sequence is valid or raise `InvalidSequenceError`.
    """
    illegal_sites = RestrictionBatch(BIOBRICK_ILLEGAL_RESTRICTION_SITES)
    results = illegal_sites.search(sequence)

    for restriction_site, match_positions in results.items():

        if len(match_positions) > 0:
            raise InvalidBioBrickError(
                '{} site found in sequence.'.format(restriction_site))

    return True


class BioBrick3AAssembler(AssemblerBase):
    """Assemble BioBrick standard parts using the 3A (three antibiotic) assembly protocol.

    References:
        http://parts.igem.org/Help:Assembly/3A_Assembly
        http://parts.igem.org/Help:Protocols/3A_Assembly
    """

    def __init__(self, config_options):
        self.plasmid_backbone_name = config_options.get('plasmid_backbone', 'pSB1C3')

        self.constant_sequences = load_biobrick_constant_sequences()

    def digest_backbone(self, backbone_name):
        """Prepare a plasmid for recombination to construct new composite BioBrick part."""
        backbone_sequence = self.constant_sequences[backbone_name]

        # pylint: disable=C0103

        # Linearize the backbone by cutting at EcoRI and then excise PstI
        linear_backbone_plus_PstI = Restriction.EcoRI.catalyze(backbone_sequence, linear=False)[0]
        _, linear_backbone = self.cut_and_validate_fragments(
            Restriction.PstI,
            linear_backbone_plus_PstI)

        return linear_backbone

    def digest_part_sequence(
        self,
        left_part_sequence : Seq,
        right_part_sequence : Seq
    ) -> Tuple[Seq, Seq]:
        """Perform the restriction digest of the left and right part.

        The left part sample is cut out with EcoRI and SpeI.
        The right part sample is cut out with XbaI and PstI.

        Return a tuple of part sequences: left digested part an right digested part.
        """
        # pylint: disable=C0103


        # Perform Left Part Restriction Digest
        (_, left_partial_part_fragment) = self.cut_and_validate_fragments(
            Restriction.EcoRI,
            left_part_sequence)

        (left_digested_part_sequence, _) = self.cut_and_validate_fragments(
            Restriction.SpeI,
            left_partial_part_fragment)

        # Perform Right Part Restriction Digest
        (_, right_partial_part_fragment) = self.cut_and_validate_fragments(
            Restriction.XbaI,
            right_part_sequence)

        (right_digested_part_sequence, _) = self.cut_and_validate_fragments(
            Restriction.PstI,
            right_partial_part_fragment)

        return left_digested_part_sequence, right_digested_part_sequence

    def cut_and_validate_fragments(self, restriction_enzyme, sequence, expected_fragments = 2):
        """Cut a sequence with a particular restriction enzyme."""
        restriction_fragments = restriction_enzyme.catalyze(sequence)
        if len(restriction_fragments) != expected_fragments:
            raise InvalidBioBrickError('Unexpected restriction cut; part is not a valid biobrick')

        return restriction_fragments

    def assemble_part_tuple(self, left_part, right_part):
        """Given two parts, assemble them using the BioBrick 3A Assembly Method.

        Method (http://parts.igem.org/Help:Protocols/3A_Assembly)
        1. Restriction digest the two parts and the new backbone:
            The left part sample is cut out with EcoRI and SpeI.
            The right part sample is cut out with XbaI and PstI.

            The linearized plasmid backbone is a linear piece of DNA and should
            be cut with EcoRI and PstI.

        2. Perform Ligation Reaction
            An equimolar quantity of all 3 restriction digest products are
            combined in a ligation reaction.

            Goal: The desired result is the left part sample's SpeI overhang
            ligated with the right part sample's XbaI overhang resulting in a
            scar that cannot be cut with any of our enzymes.

            The new composite part sample is ligated into the construction
            plasmid backbone at the EcoRI and PstI sites. When the ligation is
            transformed into cells and grown on plates with antibiotic C, only
            colonies with the correct construction survive.

        """
        left_digested_part_sequence, right_digested_part_sequence = self.digest_part_sequence(
            left_part.get_sequence().seq,
            right_part.get_sequence().seq)
        digested_backbone_sequence = self.digest_backbone(self.plasmid_backbone_name)

        fragments = [
            Dseqrecord(left_digested_part_sequence),
            Dseqrecord(right_digested_part_sequence),
            Dseqrecord(digested_backbone_sequence)
        ]

        assembly = Assembly(fragments)
        print(assembly)

    def assemble(self, assembly_requests : List[AssemblyDeclaration]) -> AssemblyList:
        """
        Strategy:
            assemblies is a collection ordered list of parts

            [
                P1 ; P2 ; P3 ; P4 ; P5
                P1 ; P2 ; P6 ; P4 ; P5
                P1 ; P2 ; P7 ; P8 ; P5
            ]

        Base case:
            Given two parts, select left (rel to 3') part and make it the left part
            and make right the right part

        Complicated:
            build an agaceny graph of parts in all assemblies
            Loop until no edges left:
                Select the edge with the most number of occurrences
                remove that edge and replace with combined part

        """

        assemblies : List[Assembly] = []
        for assembly_idx, assembly_request in enumerate(assembly_requests):
            parts = assembly_request.get_parts()

            p1 = parts[0]
            p2 = parts[1]

            print('ASSEMBLE', type(p1), p1.identifier, p2.identifier)

            # TODO: NEED TO CONFIRM PART ARE VALID BIO BRICKS -
            # ie have prefix and suffix and no bad cut sites
            check_is_valid_biobrick(p1.get_sequence().seq)

            self.assemble_part_tuple(p1, p2)
            # validate that each part has the appropriate biobrick prefix/suffix
            # confirm that each part does not contain disallowed restriction sites.


class BioBrickPartProvider(PartProvider):
    """A Part provider for returning constant regions from the BioBrick standard.

    More about the BioBrick standard can be found here:
        http://parts.igem.org/Help:Standards/Assembly/RFC10

    To configure this part provider add the following to `supergsl-config.json`:
    ```
    {
        "name": "biobrick",
        "provider_class": "supergsl.plugins.igem.BioBrickPartProvider"
    }

    References
    * http://parts.igem.org/Help:Standards/Assembly/RFC10
    """

    def __init__(self, name : str, settings : dict):
        self.name = name
        self._cached_parts: Dict[str, Part] = {}


    def get_part(self, identifier : str) -> Part:
        """Retrieve a part by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """

        try:
            return self._cached_parts[identifier]
        except KeyError:
            pass

        if identifier == 'prefix':
            reference_sequence = BIOBRICK_PREFIX_SEQUENCE
            description = 'Biobrick RFC[10] prefix'
        elif identifier == 'suffix':
            reference_sequence = BIOBRICK_SUFFIX_SEQUENCE
            description = 'Biobrick RFC[10] suffix'
        else:
            raise PartNotFoundError('Part %s not found.' % identifier)

        start = SeqPosition.from_reference(
            x=0,
            rel_to=THREE_PRIME,
            approximate=False,
            reference=reference_sequence
        )
        end = start.get_relative_position(
            x=len(reference_sequence))

        part = Part(
            identifier,
            start,
            end,
            provider=self,
            description=description,
            roles=[
                SO_HOMOLOGOUS_REGION
            ])

        self._cached_parts[identifier] = part
        return part
