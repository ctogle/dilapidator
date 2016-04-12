import pstats,cProfile,os,time

__doc__ = '''some functions for measuring function calls'''

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


