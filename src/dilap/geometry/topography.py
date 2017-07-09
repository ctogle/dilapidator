from .vec3 import vec3
from .polymath import contract, ebdxy, bintbxy, binbxy
from dilap.topology import tree
from dilap.core.plotting import *
import pdb


class topography(tree):


    '''A topography is a tree where each vertex has an associated polygon
    It can identify an x,y location with a vertex by testing containment
    of the location within the polygon associated with the vertex

    The loops can be nested with respect to containment and offzet in 
    the z direction to represent a topographical map for terrain.'''


    def avert(self,par,loop,**kws):
        '''A loop must be provided as a keyword to avert'''
        v = tree.avert(self,par,loop = loop,holes = kws.get('holes',[]))
        return v


    def loc(self,p):
        '''Determine the one or two loops which locate this point topologically;
        find the vertex such that p is in the loop but none of its childrens'''
        lst = None
        unfn = [self.root]
        while unfn:
            v = unfn.pop(0)
            if p.onbxy(v.loop):

                chs = self.below(v)
                unfn.extend(chs)
                lst = v,chs

                #return v,[]

            elif p.inbxy(v.loop):
                chs = self.below(v)
                unfn.extend(chs)
                lst = v,chs
        if lst is None:
            lst = self.root,[]
        return lst


    def locloop(self,l,override = 0,mingrad = 1.0):
        '''Determine the one or two loops which locate this loop topologically
        find the vertex such that l is in the loop but none of the children
        Potentially modify loops based on intersections and override'''

        # modify the loop of v (and its children) to respect new loop l
        def correct_loops(v,l):
            lswollen = [p.cp().ztrn(v.loop[0].z-p.z) for p in l]
            #lswollen = contract(lswollen,abs(l[0].z-v.loop[0].z)/mingrad)

            ####???#import matplotlib.pyplot as plt
            #ax = plot_axes(300)
            #ax = plot_polygon(lswollen,ax,lw = 3,col = 'b')
            #ax = plot_polygon(v.loop,ax,lw = 3,col = 'r')
            #plt.show()

            newloops = ebdxy(v.loop,lswollen)
            v.loop = newloops[0]
            if len(newloops) > 1:
                for nl in newloops[1:]:
                    nv = self.avert(self.above(v),loop = nl)
            for ch in self.below(v):
                if bintbxy(lswollen,ch.loop):
                    correct_loops(ch,lswollen)

        lst = None
        unfn = [self.root]
        while unfn:
            v = unfn.pop(0)
            if binbxy(l,v.loop):
                chs = self.below(v)
                unfn.extend(chs)
                lst = v,chs
            elif bintbxy(l,v.loop):
                # modify the existing loop to accomodate the new loop
                if override == 1:
                    correct_loops(v,l)
                    chs = self.below(v)
                    unfn.extend(chs)
                    lst = self.above(v),chs
                # modify the new loop to accomodate the existing loop
                elif override == 2:
                    raise NotImplemented
                    return v,[]
                # modify neither loop (likely causes poor results)
                else:return v,[]
        if lst is None:
            lst = self.root,[]
        return lst


    def plot(self,ax = None,l = 300,s = (vec3(-10000,0,0),vec3(10000,0,0))):
        if ax is None:ax = plot_axes(l)
        def pv(v,ax):
            plot_polygon(v.loop,ax,lw = 2)
            vp = v.loop[0]
            for c in self.below(v):
                cp = c.loop[0]
                plot_points((vp,cp),ax)
                plot_edges((vp,cp),ax,lw = 3,col = 'b')
                pv(c,ax)
        pv(self.root,ax)
        return ax
        

    def plotxy(self,ax = None,l =300):
        if ax is None:ax = plot_axes_xy(l)
        def pv(v,ax):
            ax = plot_polygon_xy(v.loop,ax,lw = 2)
            for c in self.below(v):pv(c,ax)
        pv(self.root,ax)
        return ax
