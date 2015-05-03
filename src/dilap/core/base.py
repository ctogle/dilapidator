

class base(object):

    def _def(self,key,dval,**kwargs):
        if not hasattr(self,key):
            if key in kwargs.keys():aval = kwargs[key]
            else:aval = dval
            self.__dict__[key] = aval


