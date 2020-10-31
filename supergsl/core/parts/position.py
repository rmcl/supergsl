from supergsl.core.constants import THREE_PRIME, FIVE_PRIME
from supergsl.core.exception import PartSliceError

class SeqPosition(object):
    """Store relative positions in a reference sequence.

    "Hierarchical parts" are subsequences of parent parts or reference genomes.
    The `SeqPosition` class keeps track of a position relative to the reference
    sequence.
    """

    @classmethod
    def from_reference(cls, x, rel_to, approximate, reference):
        """Create a SeqPosition object from a reference sequence."""

        # If the rel_to FIVE PRIME then reverse coordinates to make it
        # relative to the three prime position.
        if rel_to == FIVE_PRIME:
            x = len(reference) - x

        return SeqPosition(
            x,
            approximate=approximate,
            reference=reference)


    def get_relative_position(self, x, approximate=False):
        """Retrieve a new `SeqPosition` relative to the current position.

        By using the `get_relative_position` method one can instantiate a new `SeqPosition`
        object that maintains a position relative to the current `SeqPosition` object.

        For example, say you have SeqPosition corresponding to the start of the HO gene
        from s. cerevisiae. To retrieve the 500 bp immediately upstream of the gene you
        could do the following:

        ```
        gHO_start = SeqPosition(...)
        pHO_start = gHO_start.get_relative_position('ThreePrime', -500, approximate=True)
        ```
        """
        return SeqPosition(
            x=x,
            approximate=approximate,
            parent=self)


    def __init__(self, x, approximate=False, reference=None, parent=None):
        """Instantiate a SeqPosition object. Don't use this constructor directly.

        Instead use `SeqPosition.from_reference` or `get_relative_position`

        """
        self.x = x
        self.approximate = approximate

        self.reference = reference
        self.parent = parent

    def __str__(self):

        return '3\'%+dbp Approx:%s (Abs Pos: %s)' % (
            self.x,
            self.approximate,
            self.get_absolute_position_in_reference()[1]
        )

    def get_absolute_position_in_reference(self):
        """Return the absolute position in the reference sequence."""
        if self.parent:
            reference, ref_relative_x = self.parent.get_absolute_position_in_reference()
            ref_relative_x += self.x

            if ref_relative_x < 0:
                raise PartSliceError(
                    'Absolute position extends beyond the 3\' end of the '
                    'reference sequence')

            if ref_relative_x > len(reference):
                raise PartSliceError(
                    'Absolute position extends beyond the 5\' end of the '
                    'reference sequence.')

            return reference, ref_relative_x

        elif self.reference:
            return self.reference, self.x

            # Thinking that if the position is on the reverse strand of the reference
            # then rather than reverse complementing here we should just pass the
            # reversed complement sequence to the SeqPosition object.
            #if self.strand == 'C':
            #    ref = ref.reverse_complement()

            # Actually lets store everything relative to THREE_PRIME
            #if self.rel_to == THREE_PRIME:
            #    return ref, self.x
            #else:
            #    return ref, len(ref) - self.x

    def check_position_compatibility(self, p2):
        """Compare to another `SeqPosition` and determine if they share the same reference sequence and strand."""

        ref1, _ = self.get_absolute_position_in_reference()
        ref2, _ = p2.get_absolute_position_in_reference()

        if ref1 != ref2:
            raise PartSliceError('Reference sequences do not match.')
