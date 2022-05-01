from typing import Optional, Dict
from inspect import getdoc


class SuperGSLType(object):
    """Base class defining types available in SuperGSL."""

    def eval(self) -> 'SuperGSLType':
        return self

    def print(self) -> str:
        """Display details about the SuperGSL object."""
        return str(self)

    @property
    def help(cls) -> Optional[str]:
        """Return helpful documentation about this type."""
        return getdoc(cls)

    def serialize(self) -> Dict:
        """Return a serializable representation of this type."""
        return {}
