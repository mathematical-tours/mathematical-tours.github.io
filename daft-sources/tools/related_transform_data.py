"""Fixed inputs for Chapter 5, matching the original illustrative experiments."""
import numpy as np
from scipy.signal import lfilter


def triangle():
    return np.r_[np.arange(8), np.arange(8, 0, -1)].astype(float)


def notch(z):
    z = np.asarray(z, dtype=complex)
    return (z*z-np.sqrt(2)*z+1)/(z*z-.9*np.sqrt(2)*z+.81)


def notch_impulse(length=128):
    impulse=np.zeros(length);impulse[0]=1
    return lfilter([1,-np.sqrt(2),1],[1,-.9*np.sqrt(2),.81],impulse)
