import dilap.core.lsystem as lsy

import dilap.core.plotting as dtl

import pdb

###############################################################################
###
###############################################################################

def apt():
    seq  = ''
    seq += 'S<0,0.5,0.5,0,0,1,0>'
    seq += 'S<0,0.5,0.5,0,1,0,0>'
    seq += 'X<0>'
    seq += 'E<1,0>'
    seq += 'E<1,2>'
    return seq

def hall():
    seq  = ''
    seq += 'S<0,0.15,0.5,0,0,1,0>'
    seq += 'S<0,0.5,0.75,0,1,0,0>'
    seq += 'S<0,0.5,0.66,0,1,0,0>'
    seq += 'S<0,0.5,0.5,0,1,0,0>'
    seq += 'X<0>'
    seq += 'E<0,1>'
    seq += 'E<1,2>'
    seq += 'E<1,3>'
    seq += 'E<1,4>'
    seq += 'I<0,'+apt()+'>'
    return seq

def building():
    seq = ''
    seq += 'L<0>'
    seq += 'S<0,0.5,0.75,0,1,0,0>'
    seq += 'S<1,0.5,0.75,0,1,0,0>'
    seq += 'I<0,'+hall()+'>'
    seq += 'I<1,'+hall()+'>'
    seq += 'V<3>'
    return seq



def myrtapt():
    seq  = ''
    seq += 'S<0,0.3,0.5,0,0,1,0>' # make living room on left
    seq += 'S<0,0.6,0.5,0,0,1,0>' # make back room on right
    seq += 'S<2,0.5,0.4,0,1,0,0>' # make kitchen on bottom
    seq += 'S<3,0.5,0.3,0,1,0,0>' # make front room on top
    seq += 'S<2,0.6,0.5,0,0,1,0>' # make bathrom on left of kitchen
    seq += 'S<5,0.5,0.7,0,1,0,0>' # make hall closest above kitchen
    seq += 'S<2,0.5,0.4,0,1,0,0>'
    seq += 'S<0,0.5,0.3,0,1,0,0>'
    seq += 'S<4,0.3,0.5,0,0,1,0>'
    seq += 'S<9,0.5,0.5,0,1,0,0>'
    seq += 'E<1,5>E<1,3>E<5,2>E<3,7>E<3,4>E<3,8>E<3,6>E<8,0>E<4,9>E<4,10>X<1>'
    seq += 'R<1,rtype,open>'
    seq += 'R<3,rtype,open>'
    seq += 'R<5,rtype,open>'
    seq += 'R<0,rtype,closed>'
    seq += 'R<10,rtype,closed>'
    return seq

def myltapt():
    seq  = ''
    seq += 'S<0,0.7,0.5,0,0,1,0>'
    seq += 'S<1,0.4,0.5,0,0,1,0>'
    seq += 'S<1,0.5,0.4,0,1,0,0>'
    seq += 'S<3,0.5,0.3,0,1,0,0>'
    seq += 'S<1,0.4,0.5,0,0,1,0>'
    seq += 'S<1,0.5,0.7,0,1,0,0>'
    seq += 'S<5,0.5,0.4,0,1,0,0>'
    seq += 'S<2,0.5,0.3,0,1,0,0>'
    seq += 'S<4,0.7,0.5,0,0,1,0>'
    seq += 'S<4,0.5,0.5,0,1,0,0>'
    seq += 'E<0,1>E<0,3>E<3,6>E<3,7>E<3,9>E<3,8>E<1,5>E<8,2>E<4,9>E<9,10>X<0>'
    seq += 'R<0,rtype,open>'
    seq += 'R<1,rtype,open>'
    seq += 'R<3,rtype,open>'
    seq += 'R<2,rtype,closed>'
    seq += 'R<10,rtype,closed>'
    return seq

def mylbapt():
    seq  = ''
    seq += 'S<0,0.7,0.5,0,0,1,0>'
    seq += 'S<1,0.4,0.5,0,0,1,0>'
    seq += 'S<1,0.5,0.6,0,1,0,0>'
    seq += 'S<1,0.5,0.7,0,1,0,0>'
    seq += 'S<3,0.4,0.5,0,0,1,0>'
    seq += 'S<3,0.5,0.3,0,1,0,0>'
    seq += 'S<5,0.5,0.6,0,1,0,0>'
    seq += 'S<2,0.5,0.7,0,1,0,0>'
    seq += 'S<1,0.7,0.5,0,0,1,0>'
    seq += 'S<1,0.5,0.5,0,1,0,0>'
    seq += 'E<0,6>E<0,4>E<7,6>E<3,4>E<4,5>E<2,4>E<1,9>E<8,2>E<4,9>E<9,10>X<0>'
    seq += 'R<0,rtype,open>'
    seq += 'R<4,rtype,open>'
    seq += 'R<6,rtype,open>'
    seq += 'R<8,rtype,closed>'
    seq += 'R<1,rtype,closed>'
    return seq

def myfloor():
    seq  = ''
    seq += 'S<0,0.5,0.5,0,1,0,0>'
    seq += 'S<0,0.4,0.5,0,0,1,0>'
    seq += 'S<0,0.33,0.5,0,0,1,0>'
    seq += 'S<1,0.5,0.5,0,0,1,0>'
    seq += 'E<3,4>E<3,2>E<3,1>E<3,0>V<3>X<3>X<3>'
    seq += 'I<0,'+myrtapt()+'>'
    seq += 'I<1,'+myrtapt()+'>'
    seq += 'I<2,'+mylbapt()+'>'
    seq += 'I<4,'+myltapt()+'>'
    return seq

def mybuilding(lvls = 3):
    seq = ''

    seq += 'X<0>'
    for lvx in range(lvls-1):
        seq += 'L<'+str(lvx)+'>'

    for lvx in range(lvls):
        seq += 'I<'+str(lvx)+','+myfloor()+'>'

    return seq

def simplebuilding(*ags,**kws):
    #seq = building()
    seq = mybuilding(*ags,**kws)
    return seq

###############################################################################
###############################################################################





