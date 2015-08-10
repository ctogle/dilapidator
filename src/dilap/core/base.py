
__doc__ = '''A base class from which all dilapidator classes inherit'''

class base(object):
    '''A base class for all other dilapidator class'''

    def _def(self,key,dval,**kwargs):
        '''If attribute "key" is not present, 
        set to dval or kwargs[key] if it exists'''
        if not hasattr(self,key):
            if key in kwargs.keys():aval = kwargs[key]
            else:aval = dval
            self.__dict__[key] = aval


