from dilap.geometry.vec3 import vec3

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import pdb

# create one model as one on an mplt axis
def build_model2(mod,**kwargs):
    ax = kwargs['ax']
    ax = dtl.plot_points(mod.pset.gps(None),ax)
    return ax



