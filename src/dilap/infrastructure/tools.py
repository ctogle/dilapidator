import dilap.core.vector as dpv
import dilap.core.tools as dpr
import dilap.core.lsystem as dls
import dilap.mesh.tools as dtl

import matplotlib.pyplot as plt
import random as rm
import pdb

def adist(a1,a2):
    da = dpr.clamp_periodic(a1-a2,0,360)
    return da if da < 180 else 360 - da

def signedadist(a1,a2):
    return a1-a2 if a1 > a2 else a2-a1

def apply_err(x,e,m):
    xpe = dpr.clamp_periodic(x+e,0,360)
    xme = dpr.clamp_periodic(x-e,0,360)
    da1,da2 = adist(xpe,m),adist(xme,m)
    if da1 > da2:return xpe
    else:return xme

def min_adist(alpha,angles,exempt = None):
    acnt = len(angles)
    if acnt < 2:return None,None
    else:
        minad = None
        minax = None
        for a2x in range(acnt):
            if a2x == exempt:continue
            ad = adist(angles[a2x],alpha)
            if minad is None or ad < minad:
                minad = ad
                minax = a2x
        return minad,minax

# given a list of angles, gradually move them to satisfy a condition
# the condition is that da for each angle be nearly equal to their avg
# da is the minimum distance to the nearest angle, which should be ~90
# ws are weights to affect the motion per nudge
def nudge(angles,ws,target = 90,error = 1):
    def acceptable(das):return abs(min(das)-target) < error
    def measure(angles):
        das = []
        dxs = []
        for x in range(acnt):
            xda,xdx = min_adist(angles[x],angles,x)
            das.append(xda)
            dxs.append(xdx)
        return das,dxs

    oas = angles[:]
    acnt = len(angles)
    if acnt == 1:return oas
    das,dxs = measure(oas)
    tries = 0
    maxtries = 100
    while not acceptable(das) and tries < maxtries:
        tries += 1
        das,dxs = measure(oas)
        tdas = [dpr.clamp(target-da,0,target) for da in das]
        es = [tda*w*0.5 for tda,w in zip(tdas,ws)]
        oas = [apply_err(oas[x],es[x],oas[dxs[x]]) for x in range(acnt)]
    if tries == maxtries:print('tries exceeded while nuding!')
    return oas










