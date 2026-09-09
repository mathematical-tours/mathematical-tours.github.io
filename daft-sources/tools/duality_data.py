"""Shared mathematical examples used in the duality chapter figures."""
import numpy as np


def walsh_matrix(n):
    return np.array([[(-1)**((i&j).bit_count()) for j in range(n)] for i in range(n)])


def chirp(n=256):
    x = np.arange(n)/n
    return x, np.cos(np.exp(4.5*x))*np.cos(np.pi*(x-.5))**2


def letter_a(n=64):
    # An analytic block letter A, independently rasterized at pixel centers.
    y, x = np.mgrid[0:n, 0:n]/n + .5/n
    center = .43*y
    legs = (np.abs(np.abs(x-.5)-center) < .085) & (y>.08) & (y<.92)
    bar = (y>.54) & (y<.66) & (np.abs(x-.5)<.34)
    cap = (y>.08) & (y<.18) & (np.abs(x-.5)<.15)
    return 1.0-(legs|bar|cap).astype(float)
