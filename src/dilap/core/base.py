import time
import appdirs
import pstats
import cProfile
import datetime
import glob
import os
import six.moves.cPickle as pickle
import pdb


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


# convenient class for saving/loading capability
class persistent_container(base):

    def save(self):
        with open(self.spath,'wb') as h:pickle.dump(self,h)
        return self


    def load(self):
        if not os.path.exists(self.spath):return False
        with open(self.spath,'rb') as h:
            scopy = pickle.load(h)
            for key in scopy.__dict__.keys():
                scopyvalue = scopy.__dict__[key]
                self.__dict__[key] = scopyvalue
        return self


    def __init__(self,savedir,name,absolute = False):
        self.savedir = savedir
        self.name = name
        if not absolute:
            sdir = os.path.join(os.getcwd(),savedir)
        else:sdir = savedir
        if not os.path.exists(sdir):os.mkdir(sdir)
        self.spath = os.path.join(sdir,name+'.pkl')


# a generator for endlessly looping through a sequence
def roundrobin(seq):
    j,l = 0,len(seq)
    while True:
        yield seq[j]
        j += 1
        if not j < l:
            j = 0


# return the path to a safe resource directory, 
# or a full path to a file therein
res_path = os.path.join(appdirs.user_data_dir(),'dilap_resources')
def resource_path(res = None):
    '''find a resource file or the directory where it should be'''
    if res is None:rpaths = [res_path[:]]
    else:
        rpaths = []
        for x in os.walk(res_path):
            for y in glob.glob(os.path.join(x[0],res)):
                rpaths.append(y)
    return rpaths


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


def nowdt():
    dt = datetime.datetime.now()
    return dt


def timestamp(outformat = '%Y-%m-%d %H:%M:%S:000',dt = None):
    if dt is None:dt = nowdt()
    ts = dt.strftime(outformat)
    return ts
