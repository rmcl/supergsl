from supergsl.core.config import settings

class PartTable(object):
    self._parts = {}

    def get(self, part_name):
        return self._parts[part_name]


class Part(object):
    def __init__(self):
        self.name = None
        self.sequence = None


class PartProvider(object):

    def get_part(self, identifier):
        """
        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """
        raise NotImplemented('Subclass to implement.')

class PartProviderFactory(object):

