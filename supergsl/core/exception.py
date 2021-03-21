"""Define custom exceptions for SuperGSL."""

class SuperGSLError(Exception):
    """Base exception for all custom SuperGSL errors."""


class ConfigurationError(SuperGSLError):
    """Raised when a plugin or other component of the system has been improperly configured."""


class ParsingError(SuperGSLError):
    """Raised when an error is encountered parsing the SuperGSL language."""


class NotFoundError(SuperGSLError):
    """Base class representing class of errors that occur when a SuperGSL object cannot be found."""


class ConfigFileNotFound(NotFoundError):
    """Attempted to load a configuration file, but it could not be found."""


class SymbolNotFoundError(NotFoundError):
    """A symbol cannot be found in the symbol table."""


class ProviderNotFoundError(NotFoundError):
    """The specific provider cannot be found.

    This most often manifests in the context of a import statement.
    """


class FunctionNotFoundError(NotFoundError):
    """The specified function cannot be found."""


class PartNotFoundError(NotFoundError):
    """Raised by `PartProvider` when a part cannot be found when imported or referenced."""


class FunctionInvokeError(SuperGSLError):
    """An Error occurred when invoking a function.

    This is often a type mismatch or a container failure.
    """


class PartSliceError(SuperGSLError):
    """An error occured during the slicing of a part."""
