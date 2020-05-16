class BackendPipelinePass(object):

    def perform(self, ast):
        raise NotImplemented('Must subclass and implement perform')
