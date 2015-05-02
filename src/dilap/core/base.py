import appdirs,numpy,os,pdb

class base(object):

    def _def(self,key,dval,**kwargs):
        if not hasattr(self,key):
            if key in kwargs.keys():aval = kwargs[key]
            else:aval = dval
            self.__dict__[key] = aval

def rad(deg):return numpy.pi*deg/180.0
def deg(rad):return 180.0*rad/numpy.pi

res_path = os.path.join(appdirs.user_data_dir(),'dilap_resources')
def resource_path(res = None):
    if res is None:rpath = res_path[:]
    else:rpath = os.path.join(res_path,res)
    return rpath


