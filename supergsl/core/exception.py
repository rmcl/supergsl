class ConfigurationError(Exception):
    pass


class ParsingError(Exception):
    pass


class PartLocatorError(Exception):
    pass


class NotFoundError(Exception):
    pass

class ProviderNotFoundError(NotFoundError):
    pass


class FunctionNotFoundError(NotFoundError):
    pass


class PartNotFoundError(NotFoundError):
    pass


class FunctionInvokeError(Exception):
    pass


class PartSliceError(Exception):
    pass

class SymbolNotFoundError(NotFoundError):
    pass
