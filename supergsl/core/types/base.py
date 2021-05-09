class SuperGSLType(object):
    """Base class defining types available in SuperGSL."""

    def eval(self, node) -> 'SuperGSLType':
        return self
