from dilap.geometry import vec3
from dilap.core.plotting import *
import numpy

class basis(vec3):

    epsilon = 0.00000001

    @classmethod
    def from_r_theta(cls, r, th):
        a = r * numpy.cos(2 * th)
        b = r * numpy.sin(2 * th)
        return cls(a, b, 0)

    @classmethod
    def from_ab(cls, a, b):
        return cls(a, b, 0)

    def evalues(self):
        d = self.mag()
        return -d, d

    def evectors(self):
        if abs(self.y) < basis.epsilon:
            if abs(self.x) < basis.epsilon:
                #major, minor = vec3(0, 0, 0), vec3(0, 0, 0)
                major, minor = None, None
            else:
                major, minor = vec3(1, 0, 0), vec3(0, 1, 0)
        else:
            e1, e2 = self.evalues()
            major = vec3(self.y, e1 - self.x, 0)
            minor = vec3(self.y, e2 - self.x, 0)
        return major, minor

    def plot(self, p, ax, which=None):
        l = 1.5
        major, minor = self.evectors()
        if (which == 'major' or not which) and major:
            major.nrm()
            dp = major.cp().uscl(l)
            plot_point_xy(p + dp, ax, col='b')
            plot_edges_xy((p, p + dp), ax, col='b', lw=2)
        if (which == 'minor' or not which) and minor:
            minor.nrm()
            dp = minor.cp().uscl(l)
            plot_point_xy(p + dp, ax, col='g')
            plot_edges_xy((p, p + dp), ax, col='g', lw=2)
        plot_point_xy(p, ax, col='r')
        return ax
