class BaseModel(object):

    def move(self, **kwargs):
        raise NotImplementedError("Must define a 'move' method on the model being called")

    def get_name(self):
        return self.__class__.__name__
    name = property(get_name, None, None)
