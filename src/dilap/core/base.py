import os,time,appdirs,pstats,cProfile

__doc__ = '''
A base class from which dilapidator classes may inherit and some useful functions
'''

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

def profile_function(func_,*args,**kwargs):
    '''profile the function "func_" which 
    receives "*args" and "**kwargs" as input'''
    cProfile.runctx('func_(*args,**kwargs)',
        globals(),locals(),'profile.prof')
    s = pstats.Stats('profile.prof')
    s.strip_dirs().sort_stats('time').print_stats()
    os.remove('profile.prof')

def measure_time(func_name,func,*args,**kwargs):
    '''crudely measure the time it takes to call function "func"'''
    st = time.time()
    ret = func(*args, **kwargs)
    en = time.time()
    took = en-st
    return ret,took

# read a sequence until proper ">" character is found
# NOTE: handles nested subsequences...
def seqread(seq,sx):
    score = 1
    sx += 1
    while score > 0:
        sx += 1
        if seq[sx] == '<':score += 1
        elif seq[sx] == '>':
            if score > 0:score -= 1
    return sx



