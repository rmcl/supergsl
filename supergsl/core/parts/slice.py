from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.constants import (
    PART_SLICE_POSTFIX_START,
    PART_SLICE_POSTFIX_END
)
from supergsl.core.ast import SlicePosition, Part as AstPartNode
from supergsl.core.parts import Part


class ResolvePartSlicePass(BreadthFirstNodeFilteredPass):

    def get_node_handlers(self):
        return {
            'Part': self.visit_part_node,
        }

    def convert_slice_position_to_seq_position(
        self,
        parent_part: Part,
        slice_position: SlicePosition
    ):
        """Convert SuperGSL `ast.SlicePosition` to `parts.SeqPosition` relative
        to the resolved parent part.

        ast.SlicePosition has the following properties:
            - index - The position on the part sequence.
            - postfix - "S" or "E" - Whether the index is relative to the start
                (S) or end (E) of the part.
            - approximate - Boolean flag determines if the index is approximate
                (True) or exact (False).
        """
        if slice_position.postfix == PART_SLICE_POSTFIX_START or slice_position.postfix is None:
            # the index is relative to the start of the part.
            return parent_part.start.get_relative_position(
                slice_position.index,
                slice_position.approximate)
        elif slice_position.postfix == PART_SLICE_POSTFIX_END:
            # the index is relative to the end of the part.
            return parent_part.end.get_relative_position(
                slice_position.index,
                slice_position.approximate)

        raise Exception('Unknown slice postfix: "%s"' % slice_position.postfix)


    def visit_part_node(self, node: AstPartNode):
        """Visit each `ast.Part` node and perform part slicing if neccessary."""

        if not node.slice:
            # This node does not require slicing. no slice has been specified.
            return node

        parent_part = node.part
        node.parent_parts = [node.part]

        start = self.convert_slice_position_to_seq_position(parent_part, node.slice.start)
        end = self.convert_slice_position_to_seq_position(parent_part, node.slice.end)

        child_identifier = '%s[%s]' % (
            node.identifier,
            node.slice.get_slice_str()
        )
        node.part = parent_part.get_child_part_by_slice(
            child_identifier, start, end)

        if node.invert:
            raise NotImplementedError('Inverted parts not implemented yet!')

        return node
