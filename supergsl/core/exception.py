class SuperGSLError(Exception):
    pass

class ConfigurationError(SuperGSLError):
    pass


class ParsingError(SuperGSLError):
    pass


class PartLocatorError(SuperGSLError):
    pass


class NotFoundError(SuperGSLError):
    pass

class ProviderNotFoundError(SuperGSLError):
    pass


class FunctionNotFoundError(SuperGSLError):
    pass


class PartNotFoundError(SuperGSLError):
    pass


class FunctionInvokeError(SuperGSLError):
    pass


class PartSliceError(SuperGSLError):
    pass

class SymbolNotFoundError(SuperGSLError):
    pass
