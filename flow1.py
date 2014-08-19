#/usr/bin/env python
from __future__ import division
import numpy as np
import pylab as pl

pi = np.pi

A = 1   ## m^3/s
B = -1   ## m^3/s

w = 2*pi    ## rad/s
ph = pi/3   ## rad

x0 = 10     ## m

d = 1       ## m


if __name__ == "__main__":

    m1 = lambda t: A * np.cos(w * t)        ## m^3/s
    m2 = lambda t: B * np.cos(w * t + ph)   ## m^3/s

    u = lambda x,t: m1(t) / (4*pi*x**2) + m2(t) / (4*pi*(x-d)**2)

    X = np.arange(0,100,0.3)

    pl.plot(X,u(X, 0))
    pl.show()

