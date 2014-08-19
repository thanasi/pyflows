## /usr/bin/env python
## potential2d.py
##
## 2d potential flow class

from __future__ import division
import numpy as np

pi = np.pi

def zero(x,y):
    if np.isscalar(x):
        if np.isscalar(y):
            return 0
        else:
            return np.zeros_like(y)
    else:
        return np.zeros_like(x)

def one(x,y):
    if np.isscalar(x):
        if np.isscalar(y):
            return 1
        else:
            return np.ones_like(y)
    else:
        return np.ones_like(x)


class FlowField(object):

    def __init__(self, xwin, ywin, res):
        self.set_win(xwin, ywin, res)
        self.reset()

    def set_win(self, xwin, ywin, res):
        self.xmin = xwin[0]
        self.xmax = xwin[1]
        self.ymin = ywin[0]
        self.ymax = ywin[1]
        try:
            self.xres = res[0]
            self.yres = res[1]
        except TypeError:
            self.xres = res
            self.yres = res

    def reset(self):
        self.X, self.Y = np.mgrid[self.xmin:self.xmax:self.xres,
                                  self.ymin:self.ymax:self.yres]

        self.U = np.zeros_like(self.X)
        self.V = np.zeros_like(self.X)
        self.PHI = np.zeros_like(self.X)
        self.PSI = np.zeros_like(self.X)

        self.flow_objects = []

    def add_object(self,f):
        self.flow_objects.append(f)
        self.U += f.u(self.X, self.Y)
        self.V += f.v(self.X, self.Y)
        self.PHI += f.phi(self.X, self.Y)
        self.PSI += f.psi(self.X, self.Y)

    def add_point(self, x0, y0, Q):
        self.add_object(FlowPoint(x0, y0, Q))

    def add_doublet_x(self, x0, y0, Q):
        self.add_object(FlowDoubletX(x0,y0,Q))

    def add_uniform(self, angle, vel):
        self.add_object(FlowUniform(angle, vel))

class FlowObject(object):
    def __init__(self, x0, y0, Q):
        self.x0 = x0
        self.y0 = y0
        self.Q = Q

    def __repr__(self):
        return "Flow Object @ (%f,%f) of strength %f" % (self.x0, self.y0, self.Q)

    def u(self, x,y):
        """ x-component of velocity """
        return zero(x,y)

    def v(self, x,y):
        """ y-component of velocity """
        return zero(x,y)

    def phi(self, x,y):
        """ velocity potential """
        return zero(x,y)

    def psi(self, x,y):
        """ stream function """
        return zero(x,y)

class FlowPoint(FlowObject):

    def __repr__(self):
        t = "source" if self.Q>0 else "sink"
        return "Point %s @ (%f, %f) of strength %f" % (t, self.x0, self.y0, self.Q)

    def u(self, x, y):
        return one(x,y) * self.Q/(2*pi) * (x-self.x0) / ((x-self.x0)**2 + (y-self.y0)**2)

    def v(self, x, y):
        return one(x,y) * self.Q/(2*pi) * (y-self.y0) / ((x-self.x0)**2 + (y-self.y0)**2)

    def phi(self, x, y):
        return one(x,y) * self.Q/(2*pi) * np.log(np.sqrt((x-self.x0)**2 + (y-self.y0)**2))

    def psi(self,x,y):
        return one(x,y) * self.Q/(2*pi) * np.arctan2(y-self.y0,x-self.x0)

class FlowDoubletX(FlowObject):

    def __repr__(self):
        return "Doublet @ (%f, %f) of strength %f" % (self.x0, self.y0, self.Q)

    def u(self, x, y):
        return one(x,y) * -self.Q/(2*pi) * ((y-self.y0)**2 - (x-self.x0)**2) / (((x-self.x0)**2 + (y-self.y0)**2))**2

    def v(self, x, y):
        return one(x,y) * self.Q/(2*pi) * 2*(y-self.y0)*(x-self.x0) / (((x-self.x0)**2 + (y-self.y0)**2))**2

    def phi(self, x, y):
        return one(x,y) * self.Q/(2*pi) * (x-self.x0) / ((x-self.x0)**2 + (y-self.y0)**2)

    def psi(self, x,y):
        return one(x,y) * -self.Q/(2*pi) * (y-self.y0) / ((x-self.x0)**2 + (y-self.y0)**2)

class FlowUniform(FlowObject):

    def __init__(self, angle, vel):
        self.angle = angle
        self.vel =vel

    def __repr__(self):
        return "Uniform flow with velocity %f at angle %f" % (self.vel, self.angle)

    def u(self, x, y):
        return one(x,y) * self.vel * np.cos(self.angle)

    def v(self, x, y):
        return one(x,y) * self.vel * np.sin(self.angle)

    def phi(self, x, y):
        return one(x,y) * self.vel * (x*np.cos(self.angle) + y*np.sin(self.angle))

    def psi(self,x,y):
        return one(x,y) * self.vel * (-x*np.sin(self.angle) + y*np.cos(self.angle))

if __name__ == "__main__":
    # import matplotlib
    # matplotlib.use('Qt4Agg')
    import matplotlib.pyplot as pl
    from matplotlib.widgets import Slider

    NN = 200
    Q = 1
    X0 = .001
    R = .25

    XWIN = [-0.25,0.25]
    YWIN = [0,0.25]
    RES = 1e-3

    pipeflow = FlowField(XWIN, YWIN, RES)
    # pipeflow.add_point(-X0,-R,Q)
    # pipeflow.add_point(-X0,R,Q)
    # pipeflow.add_point(X0,-R,-Q)
    # pipeflow.add_point(X0,R,-Q)
    #
    # pipeflow.add_point(-X0,3*R,Q)
    # pipeflow.add_point(X0,3*R,-Q)
    #
    pipeflow.add_doublet_x(0,R,Q)


    pipeflow.add_uniform(0, 10)

    fig = pl.figure()
    ax1 = pl.subplot(211)
    # pl.gca().set_aspect('equal')

    ax2 = pl.subplot(212)
    # pl.gca().set_aspect('equal')

    pl.subplots_adjust(left=0.1, bottom=0.2)
    axcont = pl.axes([0.15, 0.1, 0.65, 0.03])
    scont = Slider(axcont, 'Contours', 20, 400, valinit=NN, valfmt='%0.0f')

    def update(val):
        ncont = scont.val
        pl.sca(ax1)
        pl.text(0,1.01,"Refreshing...", size="large", style="italic",
                color = "#000000", # backgroundcolor="#fff8dc",
                transform=ax1.transAxes, ha='left', va='bottom')
        fig.canvas.draw()

        ax1.clear()
        ax2.clear()

        pl.sca(ax1)
        pl.title("$\phi$")
        pl.contour(
            pipeflow.X,pipeflow.Y,pipeflow.PHI, int(ncont))
        pl.sca(ax2)
        pl.title("$\psi$")

        pl.contour(pipeflow.X,pipeflow.Y,pipeflow.PSI, int(ncont))

        fig.canvas.draw()

    scont.on_changed(update)


    # pl.figure()
    # pl.title("velocity")
    # pl.quiver(pipeflow.X[::30,::30],pipeflow.Y[::30,::30],
    #           pipeflow.U[::30,::30],pipeflow.V[::30,::30])
    # pl.gca().set_aspect('equal')

    update(NN)
    pl.show()


