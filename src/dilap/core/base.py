import os,appdirs

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
        return self.__dict__[key]

# return the path to a safe resource directory, 
# or a full path to a file therein
res_path = os.path.join(appdirs.user_data_dir(),'dilap_resources')
def resource_path(res = None):
    '''find a resource file or the directory where it should be'''
    if res is None:rpath = res_path[:]
    else:rpath = os.path.join(res_path,res)
    return rpath



