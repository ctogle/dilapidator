import sys,os,time,pdb

import dilap.core.base as db
from dilap.geometry.vec3 import vec3
import dilap.core.context as dcx
import dilap.modeling.model as dmo
import dilap.construct as dc
import dilap.worldly.world as dwo


def bldg(p = None,q = None,s = None):
    bfp = vec3(0,0,0).sq(50,15)
    bsq = test_bseq()
    bfa = dwo.blg.blgfactory()
    cx = bfa.new(p,q,s,footprint = bfp,sequence = bsq,floorheight = 5)
    return cx
def test_bseq():
    seq  = ''
    seq += 'L<0>'
    seq += 'L<1>'
    fseq = 'S<0,0.3,0.5,0,0,1,0>' # make living room on left
    fseq += 'S<0,0.3,0.5,0,0,1,0>' # make living room on left
    fseq += 'S<0,0.3,0.5,0,1,0,0>' # make living room on left
    fseq += 'R<2,rtype,closed>'
    fseq += 'C<0,0.2,0.2,0.2,0.2>'
    fseq += 'C<1,0.2,0.2,0.2,0.2>'
    fseq += 'C<3,0.2,0.2,0.2,0.2>'
    fseq += 'C<2,0.2,2,0.2,2>'
    fseq += 'E<1,0>'
    fseq += 'E<1,2>'
    fseq += 'E<2,3>'
    fseq += 'E<2,0>'
    #fseq += 'X<3>'
    seq += 'I<0,'+fseq+'>'
    seq += 'I<1,'+fseq+'>'
    #seq += 'E<1,0>'
    #seq += 'V<2>'
    #seq += 'V<6>'
    seq += 'X<2>'
    return seq


def stage():
    #cx = cube()
    #cx = bldg()
    cx = dwo.worldfactory().new()
    return cx


if __name__ == '__main__':
    s = time.time()
    wname = 'world0'
    wdir = os.path.join(os.getcwd(),wname)
    if not os.path.exists(wdir):os.mkdir(wdir)
    if 'profile' in sys.argv:
        db.profile_function(dc.world,stage,wdir)
    else:dc.world(stage,wdir)
    print('ran %s in %f seconds' % (__file__,time.time()-s))



