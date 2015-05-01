import imp,os,appdirs

def fetch_info():
    dil_pa = os.path.join(appdirs.user_data_dir(),'dilap_resources')
    uio_pa = os.path.join(dil_pa,'user_info.py')
    uinfo_module = imp.load_source('dilap.user_info',uio_pa)
    return uinfo_module.info


