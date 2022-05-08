from supergsl.core.exception import SuperGSLError

class SuperGSLLanguageError(SuperGSLError):
    pass

class ParsingError(SuperGSLLanguageError):
    """Raised when an error is encountered parsing the SuperGSL language."""


class BackendError(SuperGSLLanguageError):
    pass
