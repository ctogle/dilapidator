from dilap.core.plotting import *
from dilap.geometry import vec3
import dilap.geometry.polymath as pym
import numpy
from PIL import Image
import cv2
from skimage import filters, segmentation
from skimage.measure import label, regionprops
from scipy.interpolate import griddata
import opensimplex

class transform:

    def __init__(self, b, zmin=0, zmax=10):
        self.xmin, self.xmax = vec3(1, 0, 0).prjps(b)
        self.ymin, self.ymax = vec3(0, 1, 0).prjps(b)
        self.zmin, self.zmax = zmin, zmax
        self.dx = self.xmax - self.xmin
        self.dy = self.ymax - self.ymin
        self.dz = self.zmax - self.zmin

    @staticmethod
    def boundary(ix, iy, iz):
        ix = max(0.0, min(1.0, ix))
        iy = max(0.0, min(1.0, iy))
        iz = max(0.0, min(1.0, iz))
        return ix, iy, iz

    def wtoi(self, wx, wy, wz):
        ix = (wx - self.xmin) / self.dx
        iy = (wy - self.ymin) / self.dy
        iz = (wz - self.zmin) / self.dz
        return self.boundary(ix, iy, iz)

    def itow(self, ix, iy, iz):
        wx = self.xmin + ix * self.dx
        wy = self.ymin + iy * self.dy
        wz = self.zmin + iz * self.dz
        return wx, wy, wz

class heightmap:

    def __init__(self, tform, *imgs):
        self.tform = tform
        self.imgs = []
        self.modes = []
        for img in imgs:
            if isinstance(img, str):
                img, mode = load_image(img), 'additive'
            if isinstance(img, tuple):
                img, mode = img
            else:
                img, mode = img, 'additive'
            self.imgs.append(img)
            self.modes.append(mode)

    # ADD INTERPOLATION???
    def lookup(self, wx, wy):
        for img, mode in zip(self.imgs, self.modes):
            ry, rx = img.shape
            ix, iy, iz = self.tform.wtoi(wx, wy, 0)
            i, j = min(rx - 1, int(rx * ix)), min(ry - 1, int(ry * iy))
            iz = img[j, i]
            yield (ix, iy, iz), mode

    def z(self, wx, wy):
        wz = 0.0
        for (ix, iy, iz), mode in self.lookup(wx, wy):
            if mode == 'additive':
                rz = 255.0
                idz = iz / rz
                wdz = idz * self.tform.dz
                wz += self.tform.zmin + wdz
            elif mode == 'multiplicative':
                wz *= iz
            else:
                print('unknown heightmap image mode: "%s"' % mode)
        return wz

    def plot(self, ax, extent=100):
        img = sum(self.imgs)
        ax.imshow(img[::-1, :], cmap=plt.get_cmap('Greys_r'),
                  extent=(-extent, extent, -extent, extent))
        return ax

    @staticmethod
    def test():
        b = vec3(0, 0, 0).pring(100, 8)
        x, y = numpy.arange(0, 100, 1), numpy.arange(0, 100, 1)
        X, Y = numpy.meshgrid(x, y)
        img = numpy.sqrt((X - 50) ** 2 + (Y - 50) ** 2)
        #img = X + Y 
        h = heightmap(transform(b, -10, 10), img)
        ax = plot_axes_xy(140, (0, 0))
        h.plot(ax)
        plt.show()

def load_image(path):
    img = Image.open(path)
    img = numpy.array(img.convert('L'))
    return img

def show_images(*imgs):
    nrows = len(imgs)
    fig, ax = plt.subplots(nrows=nrows, ncols=1, figsize=(8 * nrows, 8), dpi=200)
    if nrows == 1:
        ax.imshow(imgs[0], cmap=plt.get_cmap('Greys_r'))
    else:
        for a, i in zip(ax, imgs):
            a.imshow(i, cmap=plt.get_cmap('Greys_r'))
    plt.show()

noise_gen = opensimplex.OpenSimplex(seed=0)

def noise_filter(img, w, sx, sy=None):
    sy = sy if sy else sx
    ny, nx = img.shape
    for i in range(nx):
        for j in range(ny):
            img[j, i] += w * noise_gen.noise2d(sx * (i / nx), sy * (j / ny))
    return img

def affine_filter(img, zmin, zmax):
    m = (zmax - zmin) / (img.max() - img.min())
    b = img.min()
    return m * img + b

def normalize(img):
    img = (img - img.min()) / img.max()
    return img

def steepness(img):
    xgrad, ygrad = numpy.gradient(img)
    grad_mag = numpy.sqrt(xgrad ** 2 + ygrad ** 2)
    #flat_mask = grad_mag <= filters.threshold_otsu(grad_mag)
    return grad_mag

def perlin(rx, ry):
    img = numpy.zeros((rx, ry))
    img = noise_filter(img, 1.0, 10.0, 10.0)
    img = noise_filter(img, 0.5, 20.0, 20.0)
    img = noise_filter(img, 0.1, 40.0, 40.0)
    img = normalize(img)
    img = cv2.blur(img, (4, 4))
    img = affine_filter(img, 0, 255)
    return img

def proximal(rx, ry, b, m=None):
    m = m if m else b
    tform = transform(b, 0, 1)
    ib = vec3(rx, ry, 1).sclps([vec3(*tform.wtoi(*p)) for p in m])
    contour = numpy.array([[[int(p.y), int(p.x)]] for p in ib])
    img = numpy.zeros((rx, ry))
    for i in range(rx):
        for j in range(ry):
            d = cv2.pointPolygonTest(contour, (j, i), True)
            img[j, i] = max(0, d) ** (0.8)
            #img[j, i] = numpy.sqrt(max(0, d)) + numpy.exp(0.1 * max(0, d / 2.0))
    return img

def pixel_polygon(mask):
    pixels = mask.astype(numpy.uint8) * 255
    pixels = cv2.dilate(pixels, (8, 8), 10)
    ret, thresh = cv2.threshold(pixels, 127, 255, 0)
    im2, contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    if contours:
        vs = cv2.approxPolyDP(contours[0], 0.1, True)
        ps = [vec3(x, y, 0) for x, y in vs[:, 0]]
        #ax = plot_axes_xy(100, o=(100, 100))
        #plot_polygon_xy(ps, ax)
        #plt.show()
        return ps

def segment_regions(img):
    clean_border = segmentation.clear_border(img)
    labeled = label(clean_border)
    regions = regionprops(labeled)
    regions = sorted(regions, key=(lambda r: r.area), reverse=True)
    masks = [(labeled == r.label) for r in regions]
    fps = []
    for m, r in zip(masks, regions):
        fp = pixel_polygon(m)
        fp = vec3(1.0 / m.shape[0], 1.0 / m.shape[1], 1).sclps(fp)
        fps.append(fp)
    return fps, masks, regions
