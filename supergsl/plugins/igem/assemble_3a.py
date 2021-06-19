"""Implement the 3A (three antibiotic) assembly protocol for BioBrick Parts."""
from typing import Dict, Tuple, List
from Bio.Seq import Seq
from Bio.Restriction import Restriction

from pydna.assembly import Assembly as PyDnaAssembly
from pydna.dseqrecord import Dseqrecord

from supergsl.core.assembly import AssemblerBase
from supergsl.core.types.assembly import (
    AssemblyDeclaration,
    AssemblyResultSet,
    Assembly
)
from .exception import InvalidBioBrickError
from .utils import (
    load_biobrick_constant_sequences,
    check_is_valid_biobrick
)

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

        assembly = PyDnaAssembly(fragments)
        print(assembly)

    def assemble(self, assembly_requests : List[AssemblyDeclaration]) -> AssemblyResultSet:
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
            parts = assembly_request.get_levels_by_factor_type('Part')

            p1 = parts[0]
            p2 = parts[1]

            print('ASSEMBLE', type(p1), p1.identifier, p2.identifier)

            # TODO: NEED TO CONFIRM PART ARE VALID BIO BRICKS -
            # ie have prefix and suffix and no bad cut sites
            check_is_valid_biobrick(p1.get_sequence().seq)

            self.assemble_part_tuple(p1, p2)
            # validate that each part has the appropriate biobrick prefix/suffix
            # confirm that each part does not contain disallowed restriction sites.

        return AssemblyResultSet(assemblies)
