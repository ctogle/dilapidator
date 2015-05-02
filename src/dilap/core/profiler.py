import pstats,cProfile,os,time

def profile_function(func_):
    cProfile.runctx('func_()',globals(),locals(),'profile.prof')
    s = pstats.Stats('profile.prof')
    s.strip_dirs().sort_stats('time').print_stats()
    os.remove('profile.prof')

def measure_time(func_name,func,*args,**kwargs):
    st = time.time()
    ret = func(*args, **kwargs)
    en = time.time()
    took = en-st
    return ret,took


