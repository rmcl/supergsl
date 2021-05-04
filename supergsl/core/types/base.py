class SuperGSLType(object):
    """Base class defining types available in SuperGSL."""

    def eval(self, ast_node : 'Node') -> 'SuperGSLType':
        """Evaluate a object when it is used in a SuperGSL program."""
        return self
