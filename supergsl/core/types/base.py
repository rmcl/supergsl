class SuperGSLType(object):
    """Base class defining types available in SuperGSL."""

    def eval(self) -> 'SuperGSLType':
        return self

    def print(self) -> str:
        """Display details about the SuperGSL object."""
        return str(self)
