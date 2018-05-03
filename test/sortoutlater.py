from dilap.geometry import vec3
from dilap.geometry import planargraph
from dilap.geometry.tensor_field import tensor_field
from dilap.core.image import transform
from dilap.core.plotting import *
import numpy

def test_tensor_field():

    def radgrid():
        elements = [
            [vec3(50, 0, 0), vec3(51, 0, 0)],
            [vec3(0, 50, 0), vec3(1, 51, 0)],
            vec3(75, 75, 0),
        ]
        return tensor_field(((e, d, w) for e in elements))

    def polyline():
        l = []
        for x in numpy.arange(0, 105, 5):
            point = vec3(x, 50 + 10 * numpy.cos(x * (2 * numpy.pi) / 50), 0)
            l.append(point)
        return tensor_field([(l, d, w)])

    def gradient(b):
        x, y = numpy.arange(0, 100, 1), numpy.arange(0, 100, 1)
        X, Y = numpy.meshgrid(x, y)
        img = numpy.sqrt((X - 50) ** 2 + (Y - 50) ** 2)
        #img = X + Y 
        tform = transform(b, -10, 10)
        elements = [((img, tform), None, w)]
        return tensor_field(elements)

    d, w = 0.01, 1.0
    b = vec3(50, 50, 0).pring(40, 8)

    #tf = radgrid()
    #tf = polyline(l)
    tf = gradient(b)

    ps = [[vec3(i, j, 0) for i in range(0, 110, 5)] for j in range(0, 110, 5)]
    ps = [x for y in ps for x in y]

    ax = plot_axes_xy(50, (50, 50))
    tf.plot(ax, ps, None)
    #plot_edges_xy(l, ax, col='r', lw=4)
    for y in range(0, 110, 20):
        stream = tf.trace((vec3(0, y, 0), vec3(1, 0, 0)), 2, 50, 'major')
        ps, dps = zip(*list(stream))
        if ps:
            plot_edges_xy(ps, ax, col='m', lw=6)
            for p in ps:
                plot_point_xy(p, ax, col='m')
    for x in range(0, 110, 10):
        stream = tf.trace((vec3(x, 0, 0), vec3(0, 1, 0)), 2, 50, 'minor')
        ps, dps = zip(*list(stream))
        if ps:
            plot_edges_xy(ps, ax, col='y', lw=6)
            for p in ps:
                plot_point_xy(p, ax, col='y')
    plot_polygon_xy(b, ax, col='r', lw=6)
    plt.show()

def test_merge_ps():
    u = planargraph()
    u.fe(u.fp(vec3(-10, 0, 0), 1), u.fp(vec3(10, 0, 0), 1))
    v = planargraph()
    v.fe(v.fp(vec3(-5, -5, 0), 1), v.fp(vec3(5, 5, 0), 1))
    w = planargraph.merge_pgs(u, v)
    w.plotxy()
    plt.show()

def test_show_images():
    from dilap.core.images import *

    path = 'someimage.png'
    show_images(
        load_image(path),
        perlin(100, 100),
        proximal(100, 100, pym.splitb(vec3(0, 0, 0).pring(10, 8), 2)),
    )
