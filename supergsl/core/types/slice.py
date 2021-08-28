from typing import Dict, Any
from .base import SuperGSLType


class Position(SuperGSLType):
    def __init__(self, index: int, postfix : str, approximate : bool):
        self.index = index
        self.postfix = postfix
        self.approximate = approximate

    def serialize(self) -> Dict[str, Any]:
        return {

            'index': self.index,
            'postfix': self.postfix,
            'approximate': self.approximate
        }

    def get_slice_pos_str(self):
        return '%s%d%s' % (
            '~' if self.approximate else '',
            self.index,
            self.postfix if self.postfix else ''
        )


class Slice(SuperGSLType):
    def __init__(self, start : Position, end : Position):
        self.start = start
        self.end = end

    @classmethod
    def from_str(cls, slice_string : str) -> 'Slice':
        # Todo: Figure out how to resolve the circular dependency without putting
        # the import here.
        from supergsl.utils.slice import parse_slice_str
        return parse_slice_str(slice_string)

    def serialize(self) -> Dict[str, Any]:
        return {
            'start': self.start.serialize(),
            'end': self.end.serialize()
        }

    def get_slice_str(self):
        return '%s:%s' % (
            self.start.get_slice_pos_str(),
            self.end.get_slice_pos_str()
        )
