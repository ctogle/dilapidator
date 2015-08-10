import imp,os,appdirs

__doc__ = '''provide access to a user info file stored in a user dependent place'''

def fetch_info():
    '''load the user info file from the user dependent resources folder'''
    dil_pa = os.path.join(appdirs.user_data_dir(),'dilap_resources')
    uio_pa = os.path.join(dil_pa,'user_info.py')
    uinfo_module = imp.load_source('dilap.user_info',uio_pa)
    return uinfo_module.info


