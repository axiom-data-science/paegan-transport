class BaseModel(object):

    def move(self, **kwargs):
        raise NotImplementedError("Must define a 'move' method on the model being called")

    @property
    def name(self):
        if hasattr(self, '_name') and self._name:
            return self._name

        return self.__class__.__name__

    @name.setter
    def name(self, value):
        self._name = value

