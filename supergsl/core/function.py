from supergsl.core.backend import DepthFirstNodeFilteredPass
from supergsl.core.exception import FunctionNotFoundException

class SuperGSLFunction(object):
    """Add a callable function to SuperGSL."""

    import_path = None
    name = None

    @classmethod
    def get_name(cls):
        if cls.name:
            return cls.name

        raise NotImplementedError('sGSL function definitions must specify a "name" in "%s".' % cls)

    @classmethod
    def get_import_path(cls):
        if cls.import_path:
            return cls.import_path

        raise NotImplementedError('sGSL function definitions must specify a "import_path" in "%s".' % cls)


    @classmethod
    def get_help(cls):
        return cls.execute.__docstr__

    @classmethod
    def get_arguments(cls):
        """Return a list of expected arguments."""
        return []

    @classmethod
    def get_return_value(cls):
        """Return the expected return value of the function."""
        return None

    def execute(self, sgsl_args):
        """Called when the function is invoke in SuperGSL."""
        pass


class FunctionSymbolTable(object):

    def __init__(self):
        self._registered_functions = {}
        self._functions = {}

    def register_function(self, function_class):
        """Register a function as available for use."""
        key = (
            function_class.get_import_path(),
            function_class.get_name()
        )
        self._registered_functions[key] = function_class

    def resolve_function(self, provider_path, function_name, alias=None):
        """Resolve a function from an import statement of a superGSL program."""

        print('Attempting to resolve function: %s, %s, %s' % (provider_path, function_name, alias))

        print(self._registered_functions)

        try:
            function_class = self._registered_functions[(provider_path, function_name)]
        except KeyError:
            raise FunctionNotFoundException('"%s" does not define a function "%s".' % (
                provider_path,
                function_name
            ))

        active_alias = alias or function_name

        if active_alias in self._functions:
            raise Exception(
                'Name conflict!! Function "%s" has already been imported.' % active_alias)

        function_inst = function_class()
        self._functions[active_alias] = function_inst

    def get_function(self, identifier):
        """Retrieve an imported function."""

        try:
            return self._functions[identifier]
        except KeyError:
            raise FunctionNotFoundException('Function "%s" not defined.' % identifier)


class InvokeFunctionPass(DepthFirstNodeFilteredPass):

    def get_node_handlers(self):
        return {
            'FunctionInvocation': self.visit_function_invoke_node,
        }
