from supergsl.core.ast import Program, Assembly, Part


class SuperGSLCoreFixtures(object):

    def get_assembly_ast(self):
        parts = [
            Part('uHO'),
            Part('pADH1'),
            Part('gERG10'),
            Part('tADH1'),
            Part('dHO'),
        ]
        assembly = Assembly(parts)

        return assembly

